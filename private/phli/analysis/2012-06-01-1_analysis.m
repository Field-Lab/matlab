% Rat glaucoma control tissue, 512 array, stained for Brn3a Brn3b

piece = '2012-06-01-1';
d01 = load_data([piece '/data001'], struct('load_neurons', true, 'load_ei', true));
d01.piece.rig = 'B';
d01 = restore_stacks(d01);
keep_vars = {'piece'; 'd01'; 'keep_vars'};


%%
% d01 = init_stacks(d01);
% d01.stacks{2,1}  = build_stack([piece '/spot/stitched.tif']);
% d01.stacks{10,1} = build_stack([piece '/confocal/rt_2012-06-01-1_Brn3ab/rt_2012-06-01-1_Brn3ab_Series021_stack.tif:0']);
% d01.stacks{11,1} = build_stack([piece '/confocal/2012-06-01-1_Brn3ab/Ser024-029_Ser034-039_stack.tif:0']);
% d01.stacks{12,2} = build_stack([piece '/confocal/rt_2012-06-01-1_Brn3ab/rt_2012-06-01-1_Brn3ab_Projection_autolevel.tif']);
% save_stacks(d01);


%% Hope this was right; live images were taken in landscape and manually rotated ccw
% d01 = array_align(d01, {2,1});
% 

%%
d01.stacks{10,1} = load_slices(d01.stacks{10,1});
stack_reslicer(d01.stacks{10,1}, 'bgpoint_color', [1 1 0]);
%%
d01.stacks{10,1}.reslice_points = gui_points;
save_stacks(d01);
clear gui_points;


%%
% d01.stacks{11,1} = load_slices(d01.stacks{11,1});
% stack_reslicer(d01.stacks{11,1}, 'bgpoint_color', [1 1 0]);
%%
% d01.stacks{11,1}.reslice_points = gui_points;
% save_stacks(d01);
% clear gui_points;

%%
% resliced = getappdata(gcf, 'resliced');
% path = '/confocal/2012-06-01-1_Brn3ab/resliced1.tif';
% imwrite(resliced, [server_data_path() piece path]);
% d01.stacks{12,1} = build_stack([piece path]);
% save_stacks(d01);


%%
d01 = array_align(d01, {2,1});

%%
d01 = stack_align(d01, {12,2}, {2,1});

%%
d01 = stack_align(d01, {12,1}, {12,2});

%%
show_electrodes_over_stack(d01, {2,1})
show_electrodes_over_stack(d01, {12,2})
show_electrodes_over_stack(d01, {12,1}, 'color', 'w')

%%
save_stacks(d01);
