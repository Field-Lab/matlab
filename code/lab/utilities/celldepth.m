function d = celldepth(C)
% CELLDEPTH     Recursively determine how many levels of nested cell arrays there are
%
% 2012-07, phli
%

if ~iscell(C), d = 0; return; end
if isempty(C), d = 1; return; end
d = 1 + max(cellfun(@celldepth, C));