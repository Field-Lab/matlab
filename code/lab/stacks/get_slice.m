function slice = get_slice(stack, islice, reload)
% GET_SLICE    Get a single slice from image stack
% usage: slice = get_slice(stack, islice, reload)
%
% inputs:   stack     Image stack struct as described in Proposal.rtf
%           islice    Index of slice, defaults to 1
%           reload    Should the slice be freshly loaded even if it has already been cached?  Defaults to false.
%
% 2010-08 phli
%

if nargin < 2
    islice = 1;
end

if nargin < 3
    reload = false;
end

slices = get_slices(stack, islice, reload);
slice = slices{1};