function list_workspaces()
% LIST_WORKSPACES
% usage: list_workspaces()
%
% See also: POP_WORKSPACE, LIST_WORKSPACES, LOAD_WORKSPACE,
% CLEAR_WORKSPACE, CLEAR_WORKSPACE_STACK
%
% Original idea from:
% http://stackoverflow.com/questions/1823668/is-there-a-way-to-push-a-matla
% b-workspace-onto-a-stack
%

c = getappdata(0, 'WORKSPACE_STACK');
for i = 1:length(c)
    fprintf('\nWorkspace %d:\n', i);
    disp(c{i});
end