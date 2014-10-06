% All 18 electrodes being used in stimulation are distributed in 3 rows.
row1 = [222 271 303 335 367 418];
row2 = [190 113 81 49 17 450];
row3 = [158 115 83 51 19 482];
elloc = [row1; row2; row3];

% creating el as a 12*3 matrix (12 sets, each set with 3 el).
el = zeros(12,3);
for i = 1:6
    el(i,:) = elloc(:,i)';
end
for i = 7:12
    if i < 10
        el(i,:) = elloc(i-6,1:3);
    else
        el(i,:) = elloc(i-9,4:6);
    end
end

% creating all possible states with 3 el. (8*3 matrix)
states = pick(0:1,3,'or');

% save the information of all stimulationg pattern into a (192 = 12*8*2) structure.
for i = 1:12
    for j = 1:8
        for k = 1:2
            idx = (i-1)*8*2 + (j-1)*2 + k;
            pattern(idx).set = i; % set ID - a number between 1 and 12
            pattern(idx).el = el(i,:); % electrode set (1*3) array
            pattern(idx).pattern = j; % state ID - a number between 1 and 8
            pattern(idx).state = states(j,:); % state (1*3) array
            pattern(idx).current = 2.5*k; % current (2.5 or 5) microamps
        end
    end
end
save('c:\documents and settings\yafi\my documents\matlab\stimulationII\pattern.mat','pattern','-mat')

%% Plotting all stim electrodes
% plot_array
% row = [row1 row2 row3];
% [elx, ely] = get_electrodeLocation(row);
% for i = 1:length(row)
%     hold on
%     plot(elx(i), ely(i),'r.','MarkerSize',30);
% end
% axis off
% saveas(figure(77),'c:\documents and settings\yafi\my documents\my dropbox\stimPlan2_elpos.fig')
