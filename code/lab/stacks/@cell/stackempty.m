function bool = stackempty(stacks, varargin)
% STACK_EMPTY   Return boolean indicating if given stack is empty
% usage: bool = stackempty(stacks, stack_spec)
% e.g.:  bool = stackempty(stacks, {2,2});
%        bool = stackempty(stacks, 2, 2);
%        bool = stackempty(stacks, 'alive_montage_hires');
%
% Returns true (i.e. empty) if stacks does not have high enough dimension
% to accomodate stack_spec.
%
% 2010-10 phli
%


bool = true;
stack_index = parse_stack_index(varargin{:});

if isempty(stack_index)
    return;
end

if any(size(stacks) < cell2mat(stack_index))
    return;
end

if isempty(stacks{stack_index{:}})
    return;
end

bool = false;