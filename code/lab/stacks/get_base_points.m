function bps = get_base_points(stack, base_index)
% GET_BASE_POINTS
% usage: bps = get_base_points(stack, base_index)
%
% 2010-08 phli
%

if isfield(stack, 'base_points') && all(size(stack.base_points) >= [base_index{:}])
    bps = stack.base_points{base_index{:}};
else
    bps = [];
end