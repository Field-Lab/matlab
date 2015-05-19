function [u,n,inds] = numcellunique(cellarr)
% NUMCELLUNIQUE     From cell array of numeric arrays, get sorted unique numeric arrays
% usage: [u,n] = numcellunique(cellarr)
%
% Output U is the set of unique numeric arrays.  Output N is the indices
% into U such that cellarr = u(n).
%
% This uses some ugly and slow kludges; a MEX version using STL list::
% would be much better but not trivial to code.
%
% 2012-02 phli
%

% Kludge attack!
strings = cellfun(@mat2str, cellarr, 'UniformOutput', false);
[ustrings, ~, n] = unique(strings);
u = cellfun(@str2num, ustrings, 'UniformOutput', false);


% If desired, build up map to numeric arrays
if nargout < 3, return; end
inds = cell(size(u));
for i = 1:length(u)
    inds{i} = find(n == i);
end