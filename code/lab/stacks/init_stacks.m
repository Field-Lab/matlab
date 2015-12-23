function stacks = init_stacks(array_id)
% INIT_STACKS   Create an initial stacks array
% usage: stacks = init_stacks([array_id])
%
% Initializes with the array placeholder stack. If array_id arg is 
% included, will also try to load the array image stack.
%
% 2010-09 phli
%

stacks = {};

array_stack = struct();
array_stack.name = 'array';
array_stack.desc = 'Placeholder for array native coordinate space.';
array_stack.paths = {};
stacks = set_stack(stacks, array_stack);

if nargin > 0
    array_image_file = get_array_image(array_id);
    ai_stack = build_stack(array_image_file, 'name', 'array_image', 'basepath', @matlab_code_path);
    stacks = set_stack(stacks, ai_stack);
end