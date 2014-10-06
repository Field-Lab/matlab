function ips = get_input_points(stack, input_index)
% GET_INPUT_POINTS
% usage: ips = get_input_points(stack, input_index)
%
% 2010-08 phli
%

input_index = parse_index(input_index);
if isfield(stack, 'input_points') && all(size(stack.input_points) >= [input_index{:}])
    ips = stack.input_points{input_index{:}};
else
    ips = [];
end