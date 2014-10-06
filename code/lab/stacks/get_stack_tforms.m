function [tform, tform_inv] = get_stack_tforms(stack, varargin)
% GET_STACK_TFORMS
% usage: [tform, tform_inv] = get_stack_tforms(stack, index)
%
% 2010-08 phli
%

index = parse_stack_index(varargin{:});

if isfield(stack, 'tforms') && all(size(stack.tforms) >= [index{:}])
    tform = stack.tforms{index{:}};
else
    tform = [];
end

if isfield(stack, 'tforms_inv') && all(size(stack.tforms_inv) >= [index{:}])
    tform_inv = stack.tforms_inv{index{:}};
else
    tform_inv = [];
end