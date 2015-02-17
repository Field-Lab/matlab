function push_workspace()
% PUSH_WORKSPACE    Save the current workspace away in memory
% usage: push_workspace()
%
% This is extremely useful if you are in the middle of working on one thing
% and need to switch gears temporarily and you don't want to get your
% workspace all gummed up and you don't want to save the current workspace
% to disk.
%
% See also: POP_WORKSPACE, LIST_WORKSPACES, LOAD_WORKSPACE,
% CLEAR_WORKSPACE, CLEAR_WORKSPACE_STACK
%
% 2011-06 phli, taken directly from:
% http://stackoverflow.com/questions/1823668/is-there-a-way-to-push-a-matlab-workspace-onto-a-stack
%

c = getappdata(0, 'WORKSPACE_STACK');
if isempty(c)
    c = {};
end

% Grab workspace
w = evalin('caller', 'whos');
names = {w.name};
s = struct;
for i = 1:numel(w)
    s.(names{i}) = evalin('caller', names{i});
end

% Push it on the stack
c{end+1} = s;
setappdata(0, 'WORKSPACE_STACK', c);