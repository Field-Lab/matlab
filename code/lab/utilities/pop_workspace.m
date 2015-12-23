function pop_workspace(num)
% POP_WORKSPACE
% usage: pop_workspace(num)
%
% See also: POP_WORKSPACE, LIST_WORKSPACES, LOAD_WORKSPACE,
% CLEAR_WORKSPACE, CLEAR_WORKSPACE_STACK
%
% 2011-05, phli, taken directly from:
% http://stackoverflow.com/questions/1823668/is-there-a-way-to-push-a-matlab-workspace-onto-a-stack
%

% Pop last workspace off stack
c = getappdata(0, 'WORKSPACE_STACK');
if isempty(c)
    warning('Nothing on workspace stack');
    return;
end

if nargin < 1
    num = length(c);
end
s = c{num};
c(num) = [];

setappdata(0, 'WORKSPACE_STACK', c);

% Do this if you want a blank slate for your workspace
evalin('caller', 'clear');

% Stick vars back in caller's workspace
names = fieldnames(s);
for i = 1:numel(names)
    assignin('caller', names{i}, s.(names{i}));
end