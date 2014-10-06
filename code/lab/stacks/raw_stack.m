function raw = raw_stack(stack)
% RAW_STACK
% usage: raw = raw_stack(stack)
%
% 2010-12 phli
%

if ~isstruct(stack)
    raw = stack;
    return
end

stack = load_slices(stack);

siz = stack_size(stack);
type = stack_datatype(stack);
raw = zeros(siz, type);

for i = 1:stack_length(stack)
    raw(:,:,:,i) = get_slice(stack, i);
end