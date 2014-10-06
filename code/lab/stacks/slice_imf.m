function imf = slice_imf(stack, islice, reload)
% SLICE_IMF    Image info IMF struct for given slice
% usage: imf = slice_imf(stack, [islice, reload])
%
% inputs:   stack   Image stack struct as described in Proposal.rtf
%           islice  Index of slice from stack, defaults 1
%           reload  Whether to reload if already cached, defaults false
%
% outputs: returns only the IMF struct for the given slice
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
    imf = fullimf(stack.pages{islice});
else
    imf = [];
end