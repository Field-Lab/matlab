%% Load datarun
piece = '2012-04-13-5';
loadopts = struct('load_params', true, 'load_neurons', true, 'load_ei', true);
d00 = load_data([piece '/data000'], loadopts);
d00.piece.rig = 'A';


%% Load existing stacks
% d00 = restore_stacks(d00, 'file', 'example_stacks.mat');


%% View stacks
imshow(get_slice(d00.stacks{...}));


%% Initialize and build new stacks
d00 = init_stacks(d00);
d00.stacks{2,1}  = build_stack([piece '/spot/6.5x/stitched.tif']);
d00.stacks{2,2}  = build_stack([piece '/spot/15x/stitched.tif']);
d00.stacks{2,3}  = build_stack([piece '/spot/6.5x/stitched_horzflip.tif']);
d00.stacks{2,4}  = build_stack([piece '/spot/15x/stitched_horzflip.tif']);
d00.stacks{10,1} = build_stack([piece '/confocal/2012-04-13-5_AnkG_TUJ1/2012-04-13-5_AnkG_TUJ1_Series027_Merged_stack_DAPI_restitched.tif:[1-9]']);
d00.stacks{10,2} = build_stack([piece '/confocal/2012-04-13-5_AnkG_TUJ1/2012-04-13-5_AnkG_TUJ1_Series064_Merged_stack_TUJ1_displaced_restitched.tif:[1-22]']);
d00.stacks{10,3} = build_stack([piece '/confocal/2012-04-13-5_DAPI/2012-04-13-5_DAPI_Series041_Merged_stack.tif:[1-14]']);
d00.stacks{10,5} = build_stack([piece '/confocal/2012-04-13-5_DAPI2/2012-04-13-5_DAPI2_Series005_stack_180rot.tif:[1-8]']);
d00.stacks{10,6} = build_stack([piece '/confocal/DAPI_stitched.tif']);
d00.stacks{11,1} = build_stack([piece '/confocal/2012-04-13-5_AnkG_TUJ1/2012-04-13-5_AnkG_TUJ1_Series080_092_Merged_stack.tif:[1-25]']);
d00.stacks{11,2} = build_stack([piece '/confocal/2012-04-13-5_AnkG_TUJ1/2012-04-13-5_AnkG_TUJ1_Series080_092_Merged_stack_enhanced_restitched.tif:[1-23]']);
d00.stacks{11,3} = build_stack(fullfile(piece, '/confocal/2012-04-13-5_AnkG_TUJ1_2/2012-04-13-5_AnkG_TUJ1_2_Series026_030_36_38_44_46_54_56_stack_merge.tif:[1-23]'));
d00.stacks{11,4} = build_stack(fullfile(piece, '/confocal/2012-04-13-5_AnkG_TUJ1_2/2012-04-13-5_AnkG_TUJ1_2_Series065_069_074_075_083_084_090_091_099_101_106_107_stack_merge.tif:[1-28]'));
d00.stacks{12,1} = build_stack([piece '/confocal/2012-04-13-5_AnkG_TUJ1/2012-04-13-5_AnkG_TUJ1_Series080_092_Merged_stac_enhanced_restitched_resliced.tif']);
d00.stacks{12,2} = build_stack(fullfile(piece, '/confocal/Tuj1_AnkG_DAPI_resliced_stitched.tif'));


%% Save stacks
% save_stacks(d00, 'file', 'example_stacks.mat');


%% Reslicing
stack = {11,2};
d00.stacks{stack{:}} = load_slices(d00.stacks{stack{:}});
stack_reslicer(d00.stacks{stack{:}}, 'points_color', 'b', 'bgpoint_color', 'g', 'method', 'linear');


%% Build alignments and tforms

%Control points, image to array coordinates
d00 = array_align(d00, {2,1});
d00 = array_align(d00, {2,2});

% Control points, image to image
d00 = stack_align(d00, {2,4}, {2,2});
d00 = stack_align(d00, {10,6}, {2,4});


% Build transform based on direct control points (image to image or image to
% array)
d00.stacks = make_stack_tform(d00.stacks, );

% Find indirect/composite string of transforms between two coordinate
% systems
d00.stacks = build_tform_routes(d00.stacks);

% Build the indirect composite transform based on routes
d00.stacks = build_tform_from_routes(d00.stacks, );

d00 = clear_tform(d00);

%% Check alignment
show_electrodes_over_stack(d00, {2,4}, 'color', 'y');
show_electrodes_over_stack(d00, {10,6}, 'color', 'y');
show_electrodes_over_stack(d00, {12,2}, 'color', 'w');

%%

d00 = calc_csd(d00);
