function slice = get_slice(stacks, stack_index, islice, reload)
% GET_SLICE    Get a single slice from image stack
% usage: slice = get_slice(stacks, stack_index, islice, reload)
%
% inputs:   stacks      Stacks cell array
%           stack_index Stack index, see 
%           islice      Index of slice, defaults to 1
%           reload      Should the slice be freshly loaded even if it has already been cached?  Defaults to false.
%
% 2010-08 phli
%

if nargin < 3
    islice = 1;
end

if nargin < 4
    reload = false;
end

index = parse_stack_index(stack_index);
stack = stacks{index{:}};
slice = get_slice(stack, islice, reload);