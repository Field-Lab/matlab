function leave(varargin)
% LEAVE    Clear all variables from workspace /except/ those given as args
%
% Usage: leave importantdata
%
% By default, will check whether all the variables given actually exist; if
% a variable doesn't exist, maybe you mistyped it, so we abort without
% clearing.  Passing '-f' as first or last arg will override this.
%
% 2010-01 phli

leave = varargin;
force = false;
if strcmp(varargin{1}, '-f')
    force = true;
    leave = leave(2:end);
end
if strcmp(varargin{end}, '-f')
    force = true;
    leave = leave(1:(end-1));
end

vars = evalin('base', 'who');  % Cell array of current workspace vars

if ~force
    % Check that leaves exist; if not there's probably a typo
    not_found = setdiff(leave, vars);
    if ~isempty(not_found)
        error('LEAVE:VariablesNotFound', ['Some variables were not found in the workspace: ' join(not_found, ', ')]);
    end
end

% Clear anything not on the leave list
clears = setdiff(vars, leave);
evalin('base', ['clear ' sprintf('%s ', clears{:})]);