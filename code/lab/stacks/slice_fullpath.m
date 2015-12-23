function fullpath = slice_fullpath(stack, islice)
% SLICE_FULLPATH    Get the full path to the given slice
% usage: fullpath = slice_fullpath(stack, islice)
%
% inputs:   stack   Image stack struct as described in Proposal.rtf
%           islice  Index of slice from stack, defaults 1
%
% 2010-08 phli
%

if nargin < 2
    islice = 1;
end

path = stack.paths{islice};
fullpath = add_base_path(path, get_stack_basepath(stack));