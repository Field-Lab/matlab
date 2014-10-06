function fullimf = slice_fullimf(stack, islice, reload)
% SLICE_FULLIMF    Get full image info IMF struct for given slice
% usage: fullimf = slice_fullimf(stack, [islice, reload])
%
% inputs:   stack   Image stack struct as described in Proposal.rtf
%           islice  Index of slice from stack, defaults 1
%           reload  Whether to reload if already cached, defaults false
%
% output: The entire IMF database for the entire multipage file that slice
% is taken from.  So if slice is from a file with 4000 pages (i.e. 4000
% stacked images) then the output will be a 4000 element struct vector.
%
% 2010-08 phli
%

if nargin < 2
    islice = 1;
end

if nargin < 3
    reload = false;
end


path = stack.paths{islice};
if ~reload && isfield(stack, 'imfs') && stack.imfs.isKey(path)
    fullimf = stack.imfs(path);
else
    fullimf = [];
end