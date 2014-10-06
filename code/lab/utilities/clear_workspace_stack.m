function clear_workspace_stack()
% CLEAR_WORKSPACE_STACK
% usage: clear_workspace_stack()
%
% See also: POP_WORKSPACE, LIST_WORKSPACES, LOAD_WORKSPACE,
% CLEAR_WORKSPACE, CLEAR_WORKSPACE_STACK
%
% Original idea from:
% http://stackoverflow.com/questions/1823668/is-there-a-way-to-push-a-matla
% b-workspace-onto-a-stack
%

rmappdata(0, 'WORKSPACE_STACK');