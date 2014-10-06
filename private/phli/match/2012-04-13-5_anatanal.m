piece = '2012-04-13-5';
rig = 'A';
loadopts = struct('load_params', true, 'load_neurons', true, 'load_ei', true);
keep_vars = {'piece'; 'loadopts'; 'keep_vars'; 'd00'};

%%
d00 = load_data([piece '/data000'], loadopts);
d00 = calc_csd(d00);
d00 = restore_stacks(d00);
d00.piece.rig = rig;

%%
% d00.piece.array_id = 1504;
% d00 = init_stacks(d00);
% d00.stacks{2,1}  = build_stack([piece '/spot/6.5x/stitched.tif'], 'basepath', '~/Desktop');
% d00.stacks{2,2}  = build_stack([piece '/spot/15x/stitched.tif'],  'basepath', '~/Desktop');
% d00.stacks{2,3}  = build_stack([piece '/spot/6.5x/stitched_horzflip.tif']);
% d00.stacks{2,4}  = build_stack([piece '/spot/15x/stitched_horzflip.tif']);
% d00.stacks{10,1} = build_stack([piece '/confocal/2012-04-13-5_AnkG_TUJ1/2012-04-13-5_AnkG_TUJ1_Series027_Merged_stack_DAPI_restitched.tif:[1-9]']);
% d00.stacks{10,2} = build_stack([piece '/confocal/2012-04-13-5_AnkG_TUJ1/2012-04-13-5_AnkG_TUJ1_Series064_Merged_stack_TUJ1_displaced_restitched.tif:[1-22]']);
% d00.stacks{10,3} = build_stack([piece '/confocal/2012-04-13-5_DAPI/2012-04-13-5_DAPI_Series041_Merged_stack.tif:[1-14]']);
% d00.stacks{10,5} = build_stack([piece '/confocal/2012-04-13-5_DAPI2/2012-04-13-5_DAPI2_Series005_stack_180rot.tif:[1-8]']);
% d00.stacks{10,6} = build_stack([piece '/confocal/DAPI_stitched.tif']);
% d00.stacks{11,1} = build_stack([piece '/confocal/2012-04-13-5_AnkG_TUJ1/2012-04-13-5_AnkG_TUJ1_Series080_092_Merged_stack.tif:[1-25]'], 'basepath', '~/Desktop');
% d00.stacks{11,2} = build_stack([piece '/confocal/2012-04-13-5_AnkG_TUJ1/2012-04-13-5_AnkG_TUJ1_Series080_092_Merged_stack_enhanced_restitched.tif:[1-23]']);
% d00.stacks{12,1} = build_stack([piece '/confocal/2012-04-13-5_AnkG_TUJ1/2012-04-13-5_AnkG_TUJ1_Series080_092_Merged_stac_enhanced_restitched_resliced.tif']);
% d00.stacks = fix_basepaths(d00.stacks);
% d00.stacks{11,3} = build_stack(fullfile(piece, '/confocal/2012-04-13-5_AnkG_TUJ1_2/2012-04-13-5_AnkG_TUJ1_2_Series026_030_36_38_44_46_54_56_stack_merge.tif:[1-23]'));
% d00.stacks{11,4} = build_stack(fullfile(piece, '/confocal/2012-04-13-5_AnkG_TUJ1_2/2012-04-13-5_AnkG_TUJ1_2_Series065_069_074_075_083_084_090_091_099_101_106_107_stack_merge.tif:[1-28]'));
% d00.stacks{12,2} = build_stack(fullfile(piece, '/confocal/Tuj1_AnkG_DAPI_resliced_stitched.tif'));
% d00.stacks{11,5} = build_stack(fullfile(piece, '/confocal/Tuj1_2_3_stitch.tif:[1-39]'));
% d00.stacks{11,3} = []; % Deprecate this in favor of 11,5
% d00.stacks{11,6} = build_stack(fullfile(piece, '/confocal/2012-04-13-5_Tuj1_AnkG_3/2012-04-13-5_Tuj1_AnkG_3_042_044_049_051_060_061_065_066_stitch.tif:[1-34]'));

save_stacks(d00);

%% Copy over reslice points to new, expanded stack
% rps = d00.stacks{11,3}.reslice_points;
% size(d00.stacks{11,3}.data{1})
% size(d00.stacks{11,5}.data{1})
% diff([size(d00.stacks{11,3}.data{1}); size(d00.stacks{11,5}.data{1})])
% offset = [0 1166];
% d00.stacks{11,5}.reslice_points = cellfun(@(p) p+zerodimtoempty(repmat(offset, size(p,1), 1)), rps, 'UniformOutput', false);
% stack_reslicer(d00.stacks{11,5}, 'points_color', 'b', 'bgpoint_color', 'g', 'method', 'linear');


%% Reslicing
stack = {11,6};
d00.stacks{stack{:}} = load_slices(d00.stacks{stack{:}});
stack_reslicer(d00.stacks{stack{:}}, 'points_color', 'r', 'bgpoint_color', 'g', 'method', 'linear');
d00.stacks{stack{:}}.reslice_points = gui_points;
save_stacks(d00);

impath = '/confocal/2012-04-13-5_Tuj1_AnkG_3/2012-04-13-5_Tuj1_AnkG_3_042_044_049_051_060_061_065_066_resliced.tif';
imwrite(gui_resliced, fullfile(server_path(), piece, impath));


%% Alignment
% d00 = array_align(d00, {2,1});
% d00 = array_align(d00, {2,2});

% Used stack align to align {2,3} to {2,1}, but actually they are just
% horizontal flips of each other so should probably set transform matrix
% manually...

% d00 = stack_align(d00, {10,6}, {2,3});
% Looks pretty good except around top and top-right

% d00 = stack_align(d00, {12,1}, {10,6});

%% Generate composite transforms as necessary
% d00.stacks = build_tform_from_routes(d00.stacks, {2,3}, {1,1});
% d00.stacks = build_tform_from_routes(d00.stacks, {10,6}, {1,1});
% d00.stacks = build_tform_from_routes(d00.stacks, {12,1}, {1,1});

%% Check electrode positions
show_electrodes_over_stack(d00, {2,4}, 'color', 'y'); % Lower res; should do alignment for 2,4
show_electrodes_over_stack(d00, {10,6}, 'color', 'y');
show_electrodes_over_stack(d00, {12,2}, 'color', 'w');
% Pretty good, although in the branches of the Y the electrodes are a
% little low in 12,1 compared to 2,3 / 2,2.


%% Mark anatomical cell centers
anatmarkfile = fullfile(server_path(), piece, 'anatmarks');
load(anatmarkfile);

stack_point_picker(d00.stacks{12,2}, 'points_color', 'g', 'bgpoint_color', 'r', 'channels', 3, 'points', anatmarks1)
anatmarks1 = anatmarks1{1};
save(anatmarkfile, '-append', 'anatmarks1');


%% GUI mark EI soma and initial segment
phlmarks = fullfile(server_path(), piece, 'data000', 'phlmarks');
load(phlmarks);

celltype = 1;
cellids = d00.cell_types{celltype}.cell_ids;
numcells = length(cellids);
load(fullfile(server_path, piece, 'data000/marks1'));

for i = 1:numcells
    cellid = cellids(i);
    [~, marks1{i}] = mark_ei_somseg(d00, cellid);
end
backup = marks1;

save(phlmarks, '-append', 'marks1');


%% Plot marks over anatomy!
ai = load_array_info(d00, 2);
stack = {12,2};

d00.stacks{stack{:}}.basepath = server_data_path();
array_tform_inv = d00.stacks{stack{:}}.tforms_inv{1,1};
resliced = get_slice(d00.stacks, {stack{:}});
resliced = resliced(:,:,1:3);
imshow(resliced); hold on

load(fullfile(server_path, piece, 'data000/marks1'));
for i = 1:length(marks1)
    mark = marks1{i};
    mark = tforminv(ai.T_array_to_vision_ei, mark);
    arrow = tformfwd(array_tform_inv, mark);

    h = imarrow(gca, arrow);
    c = get(h, 'Children');
    set(c, 'Color', 'w');
    set(c(4), 'Visible', 'off'); % Kill drop shadow
    
    mid = mean(arrow);
    h = text(mid(1), mid(2), sprintf('+P %d', d00.cell_types{1}.cell_ids(i)));
    set(h, 'Color', 'w');
    
    drawnow;
end

load(fullfile(server_path, piece, 'data000/marks2'));
for i = 1:length(marks2)
    mark = marks2{i};
    mark = tforminv(ai.T_array_to_vision_ei, mark);
    arrow = tformfwd(array_tform_inv, mark);

    h = imarrow(gca, arrow);
    c = get(h, 'Children');
    set(c, 'Color', 'w');
    set(c(3), 'LineStyle', '--');
    set(c(4), 'Visible', 'off'); % Kill drop shadow

    mid = mean(arrow);
    h = text(mid(1), mid(2), sprintf('-P %d', d00.cell_types{2}.cell_ids(i)));
    set(h, 'Color', 'w');
    
    drawnow;
end

load(fullfile(server_path, piece, 'data000/marks3'));
for i = 1:length(marks3)
    mark = marks3{i};
    mark = tforminv(ai.T_array_to_vision_ei, mark);
    arrow = tformfwd(array_tform_inv, mark);

    h = imarrow(gca, arrow);
    c = get(h, 'Children');
    set(c, 'Color', 'y');
    set(c(4), 'Visible', 'off'); % Kill drop shadow

    mid = mean(arrow);
    h = text(mid(1), mid(2), sprintf('+M %d', d00.cell_types{3}.cell_ids(i)));
    set(h, 'Color', 'y');
    
    drawnow;
end

d01ser = load_data([piece '/data001/serial/serial']);
d01ser = load_neurons(d01ser);
d01ser = load_txt_cell_types(d01ser, 'marksclassification');
load(fullfile(server_path, piece, 'data001/serial/marks4'));
for i = 1:length(marks4)
    mark = marks4{i};
    mark = tforminv(ai.T_array_to_vision_ei, mark);
    arrow = tformfwd(array_tform_inv, mark);

    h = imarrow(gca, arrow);
    c = get(h, 'Children');
    set(c, 'Color', 'y');
    set(c(3), 'LineStyle', '--');
    set(c(4), 'Visible', 'off'); % Kill drop shadow

    mid = mean(arrow);
    h = text(mid(1), mid(2), sprintf('-M01ser %d', d01ser.cell_types{12}.cell_ids(i)));
    set(h, 'Color', 'y');
    
    drawnow;
end

load(fullfile(server_path, piece, 'data000/marks5'));
for i = 1:length(marks5)
    mark = marks5{i};
    mark = tforminv(ai.T_array_to_vision_ei, mark);
    arrow = tformfwd(array_tform_inv, mark);
    
    h = imarrow(gca, arrow);
    c = get(h, 'Children');
    set(c, 'Color', 'w');
    set(c(3), 'LineStyle', '-.');
    set(c(4), 'Visible', 'off'); % Kill drop shadow
    
    mid = mean(arrow);
    h = text(mid(1), mid(2), sprintf('SBC %d', d00.cell_types{5}.cell_ids(i)));
    set(h, 'Color', 'w');
    
    drawnow;
end

load(fullfile(server_path, piece, 'data000/marks6'));
for i = 1:length(marks6)
    mark = marks6{i};
    mark = tforminv(ai.T_array_to_vision_ei, mark);
    arrow = tformfwd(array_tform_inv, mark);
    
    h = imarrow(gca, arrow);
    c = get(h, 'Children');
    set(c, 'Color', 'c');
    set(c(3), 'LineStyle', ':');
    set(c(3), 'LineWidth', 4);
    set(c(4), 'Visible', 'off'); % Kill drop shadow
    
    mid = mean(arrow);
    h = text(mid(1), mid(2), sprintf('-A %d', d00.cell_types{6}.cell_ids(i)));
    set(h, 'Color', 'c');
    
    drawnow;
end

load(fullfile(server_path, piece, 'data000/marks9'));
for i = 1:length(marks9)
    mark = marks9{i};
    mark = tforminv(ai.T_array_to_vision_ei, mark);
    arrow = tformfwd(array_tform_inv, mark);
    
    h = imarrow(gca, arrow);
    c = get(h, 'Children');
    set(c, 'Color', 'c');
    set(c(3), 'LineStyle', ':');
    set(c(3), 'LineWidth', 4);
    set(c(4), 'Visible', 'off'); % Kill drop shadow
    
    mid = mean(arrow);
    h = text(mid(1), mid(2), sprintf('+A %d', d00.cell_types{6}.cell_ids(i)));
    set(h, 'Color', 'c');
    
    drawnow;
end

load(fullfile(server_path, piece, 'data000/marks14'));
for i = 1:length(marks14)
    mark = marks14{i};
    mark = tforminv(ai.T_array_to_vision_ei, mark);
    arrow = tformfwd(array_tform_inv, mark);
    
    h = imarrow(gca, arrow);
    c = get(h, 'Children');
    set(c, 'Color', 'r');
    set(c(3), 'LineStyle', ':');
    set(c(3), 'LineWidth', 4);
    set(c(4), 'Visible', 'off'); % Kill drop shadow
    
    mid = mean(arrow);
    h = text(mid(1), mid(2), sprintf('+L %d', d00.cell_types{14}.cell_ids(i)));
    set(h, 'Color', 'r');
    
    drawnow;
end