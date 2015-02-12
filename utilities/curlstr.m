function strs = curlstr(pattern)
% CURLSTR    Generate series of strings based on pattern, as Unix curl
% usage: strs = curlstr(pattern)
%
% example: curlstr('asdf[1-2]_[0-6:2]')
%      ==> {'asdf1_0' 'asdf1_2' 'asdf1_4' 'asdf1_6' 'asdf2_0' 'asdf2_2' 'asdf2_4' 'asdf2_6'}
%
% 2010-07 phli
% ToDo: allow negative numbers?
% ToDo: allow escaping of '[]'?
% ToDo: allow character based iteration to iterate through a-z A-Z etc.?
%


% Find []
paramstart = strfind(pattern, '[');
paramend = strfind(pattern, ']');
if isempty(paramstart) || isempty(paramend)
    strs = {pattern};
    return;
end


% Segment the pattern string at []
pre   = pattern(1:paramstart(1));
param = pattern(paramstart(1):paramend(1));
post  = pattern(paramend(1):end);

% Strip the []
pre   = pre(1:end-1);
param = param(2:end-1);
post  = post(2:end);


strs = {};
params = split(param, ',');
for i = 1:length(params)
    
    % Parse the curl parameters...
    splitparams = split(params{i}, ':');

    startend = split(splitparams{1}, '-');
    numstart = str2double(startend{1});
    numpad = num2str(length(startend{1}));
    
    if length(startend) < 2
        numend = numstart;
    else
        numend = str2double(startend{2});
    end
    
    if length(splitparams) < 2
        numstep = 1;
    else
        numstep = str2double(splitparams{2});
    end    
    


    % Build up the series of curled strings
    sformat = ['%.' numpad 'd'];    
    for j = numstart:numstep:numend
        crl = sprintf(sformat, j);
        jpre = [pre crl];
        
        if isempty(post)
            strs{end+1} = jpre;
        else
            % Recursively do any additional curls after this [] and string
            % concat them to each of the current results
            posts = curlstr(post);
            
            for k = 1:length(posts)
                strs{end+1} = [jpre posts{k}];
            end
        end
    end
    
    
    
    % Hack to deal with empty brackets
    if isempty(j)        
        if isempty(post)
            strs{end+1} = pre;
        else
            % Recursively do any additional curls after this []
            posts = curlstr(post);
            
            for k = 1:length(posts)
                strs{end+1} = [pre posts{k}];
            end
        end
    end
end