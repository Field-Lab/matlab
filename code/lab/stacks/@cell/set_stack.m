function stacks = set_stack(stacks, stack, varargin)
% SET_STACK     Assign stack to stacks array
% usage: stacks = set_stack(stacks, stack, stack_spec)
%        stacks = set_stack(stacks, stack)
%
% Stack spec is as described in PARSE_STACK_INDEX.  If the stack has a
% "name" field that matches a named stack index in STACK_LABELS, then that
% value can be used and the stack_spec argument can be left off.
%
% See also: PARSE_STACK_INDEX, STACK_LABELS
%
% 2010-09 phli
%

if ~isempty(varargin)
    index = parse_stack_index(varargin{:});
elseif isfield(stack, 'name')
    index = parse_stack_index(stack.name);
else
    index = {};
end


if isempty(index)
    warning('Index not found');
    return
end


stacks{index{:}} = stack;