% This function takes a specified contrast (eg [0.48 0.48 0.48]) and a
% specified direction ('left' or 'right') and finds which stimulus code
% that corresponds to

function config_num = findStimType(datarun, contrast, direction)


if strcmp(direction , 'left')
    dir = -1;
else
    dir = 1;
end

patterns = datarun{1,2}.stimulus.trials;
RGB = zeros(length(patterns), 3);
X_DELTA = zeros(length(patterns),1);
try
for i =1:length(patterns)
    RGB(i, :) = patterns(i).RGB;
    X_DELTA(i, :) = patterns(i).X_DELTA;
end
catch
    disp('problem with RGB or X_DELTA')
end

X_DELTA  = X_DELTA./abs(X_DELTA); % make +1 for right or -1 for left

RGB_match = ismember(RGB,contrast,'rows');
X_DELTA_match = ismember(X_DELTA,dir,'rows');

match = [RGB_match, X_DELTA_match];
stim_matches = ismember(match,[1,1],'rows');
matches = find(stim_matches == 1);
config_num = datarun{1,2}.stimulus.trial_list(matches(1));

end
