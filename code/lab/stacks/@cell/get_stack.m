function stack = get_stack(stacks, varargin)
% GET_STACK     Flexible indexing into stacks cell array
% usage: stack = get_stack(stacks, {2,1})
%        stack = get_stack(stacks, [2,1])
%        stack = get_stack(stacks, 2, 1)
%        stack = get_stack(stacks, 'alive_montage_lores')
%
% 2010-08 phli
%

index = parse_stack_index(varargin{:});

if isempty(index)
    stack = struct([]);
end

stack = stacks{index{:}};