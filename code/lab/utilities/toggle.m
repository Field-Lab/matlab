function opposite = toggle(val)
% TOGGLE    Return the opposite value
% usage: opposite = toggle(val)
%
% Works on logicals, and on the strings 'on'/'off' and 'true'/'false'.
% Works on doubles 0 1, although logicals are preferred.
%
% 2011-07 phli

opposite = [];
switch class(val)
    case 'logical'
        opposite = ~val;
    case 'char'
        if strcmp(val, 'on')
            opposite = 'off';
        elseif strcmp(val, 'off')
            opposite = 'on';
        elseif strcmp(val, 'true')
            opposite = 'false';
        elseif strcmp(val, 'false')
            opposite = 'true';
        end
    case 'double'
        if val == 0
            opposite = 1;
        elseif val == 1
            opposite = 0;
        end
end