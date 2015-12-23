%% Stacks
stack  = build_stack('2011-02-02/HCN1/HCN1_Series025_Merged_stack.tif:[1-10]');
stack.basepath  = '/snle/data';
imshow(get_slice(stack));
stack = load_slices(stack);
stack_reslicer(stack, 'points', gui_points, 'points_color', [1 0 1]);
resliced = getappdata(gcf, 'resliced');