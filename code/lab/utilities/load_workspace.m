function load_workspace(ind)
% LOAD_WORKSPACE
% usage: load_workspace(ind)
%
% See also: POP_WORKSPACE, LIST_WORKSPACES, LOAD_WORKSPACE,
% CLEAR_WORKSPACE, CLEAR_WORKSPACE_STACK
%
% Original idea from:
% http://stackoverflow.com/questions/1823668/is-there-a-way-to-push-a-matla
% b-workspace-onto-a-stack
%

% Get workspace from stack
c = getappdata(0, 'WORKSPACE_STACK');
if ind > length(c)
    warning('Workspace not found');
    return;
end
s = c{ind};

% Stick vars back in caller's workspace
names = fieldnames(s);
for i = 1:numel(names)
    assignin('caller', names{i}, s.(names{i}));
end