function clear_workspace(ind)
% CLEAR_WORKSPACE
% usage: clear_workspace(ind)
%
% See also: POP_WORKSPACE, LIST_WORKSPACES, LOAD_WORKSPACE,
% CLEAR_WORKSPACE, CLEAR_WORKSPACE_STACK
%
% Original idea from:
% http://stackoverflow.com/questions/1823668/is-there-a-way-to-push-a-matla
% b-workspace-onto-a-stack
%

c = getappdata(0, 'WORKSPACE_STACK');
if ind > length(c)
    warning('Workspace not found');
    return;
end
c(ind) = [];
setappdata(0, 'WORKSPACE_STACK', c);