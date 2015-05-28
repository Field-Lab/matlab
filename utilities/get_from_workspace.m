function get_from_workspace(varargin)
% GET_FROM_WORKSPACE
% usage1: get_from_workspace var1 var2 var3
% usage2: get_from_workspace(workspace_nums, var1, var2, var3)
%
% 2011-07 phli
%

c = getappdata(0, 'WORKSPACE_STACK');
if isempty(c)
    warning('Nothing on workspace stack');
    return;
end


if isnumeric(varargin{1})
    % usage2
    workspace_nums = varargin{1};
    vars = varargin(2:end);
else
    % usage1
    workspace_nums = 1:length(c);
    vars = varargin;
end

for ind = fliplr(workspace_nums)
    if ind > length(c)
        warning(['Workspace ' num2str(ind) ' not found.']);
    end
    
    s = c{ind};
    for j = 1:length(vars)
        varname = vars{j};
        if isfield(s, varname)
            assignin('caller', varname, s.(varname));
        end
    end
end