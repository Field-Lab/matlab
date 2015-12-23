piece = '2013-05-28-5';
rig = 'A';
loadopts = struct('load_params', true, 'load_neurons', true, 'load_ei', true);
keep_vars = {'piece'; 'loadopts'; 'keep_vars'; 'd00'};


%%
d00cm = load_data([piece '/d00-01-02-r/data000/data000'], loadopts);
d00cm = calc_csd(d00cm);
d00cm = restore_stacks(d00cm);
d00cm.piece.rig = rig;


%%
% d00cm.piece.array_id = 1504;
% d00cm = init_stacks(d00cm);
% d00cm.stacks{2,1}  = build_stack([piece '/Spot/4x/stitched.tif']);
% d00cm.stacks{2,2}  = build_stack([piece '/Spot/15x/stitched.tif']);
% d00cm.stacks{2,3}  = build_stack([piece '/Spot/15x/stitched_90CCW.tif']);
% d00cm.stacks{10,1}  = build_stack([piece '/confocal/good images/20x/DAPI-AF647_stack.tif:[1-25]']);
% d00cm.stacks{10,2}  = build_stack([piece '/confocal/2013-07-19/20x/DAPI_lower_area_stack.tif:[1-10]']);
% d00cm.stacks{10,3}  = build_stack([piece '/confocal/2013-07-31/tile20x_7_top_stack.tif:[1-8]']);
% d00cm.stacks{10,4}  = build_stack([piece '/confocal/2013-07-31/tile20x_7_bottom_stack.tif:[1-8]']);
% d00cm.stacks{10,5}  = build_stack([piece '/confocal/2013-07-31/tile20x_8_stack.tif:[1-12]']);
% d00cm.stacks{10,6}  = build_stack([piece '/confocal/2013-07-31/tile20x_9_top_stack.tif:[1-8]']);
% d00cm.stacks{10,7}  = build_stack([piece '/confocal/2013-07-31/tile20x_9_bottom_stack.tif:[1-9]']);
% d00cm.stacks{11,1}  = build_stack([piece '/confocal/good images/20x/DAPI_-AF647_stack_resliced.tif']);
% d00cm.stacks{11,2}  = build_stack([piece '/confocal/2013-07-19/20x/DAPI_lower_area_stack_resliced.tif']);
% d00cm.stacks{11,3}  = build_stack([piece '/confocal/2013-07-31/tile20x_7_top_stack_resliced.tif']);
% d00cm.stacks{11,4}  = build_stack([piece '/confocal/2013-07-31/tile20x_7_bottom_stack_resliced.tif']);
% d00cm.stacks{11,5}  = build_stack([piece '/confocal/2013-07-31/tile20x_8_stack_resliced.tif']);
% d00cm.stacks{11,6}  = build_stack([piece '/confocal/2013-07-31/tile20x_9_top_stack_resliced.tif']);
% d00cm.stacks{11,7}  = build_stack([piece '/confocal/2013-07-31/tile20x_9_bottom_stack_resliced.tif']);
% d00cm.stacks{11,8}  = build_stack([piece '/confocal/DAPI_20x_stitch_resliced.tif']);
% d00cm.stacks{12,1}  = build_stack([piece '/confocal/2013-07-30/tile63x_3_4_5_merged.tif:[1-15]']);
% d00cm.stacks{12,2}  = build_stack([piece '/confocal/2013-08-02/tile63x_10_11_stack.tif:[1-15]']);
% d00cm.stacks{12,3}  = build_stack([piece '/confocal/2013-08-06/bottom_tip_merged_stack.tif:[1-16]']);
% d00cm.stacks{12,4}  = build_stack([piece '/confocal/2013-08-06/tile63x_13_14_stack.tif:[1-19]']);
% d00cm.stacks{12,5}  = build_stack([piece '/confocal/2013-08-06/tile63x_15_16_stack.tif:[1-12]']);
% d00cm.stacks{12,6}  = build_stack([piece '/confocal/2013-08-06/tile63x_17_18_stack.tif:[1-13]']);
% d00cm.stacks{13,1}  = build_stack([piece '/confocal/2013-07-30/tile63x_3_4_5_resliced.tif']);
% d00cm.stacks{13,2}  = build_stack([piece '/confocal/2013-08-02/tile63x_10_11_resliced.tif']);
save_stacks(d00cm);


%% Alignment
% d00cm = array_align(d00cm, {2,1});
% d00cm = array_align(d00cm, {2,2});
% d00cm = array_align(d00cm, {2,3});
% show_electrodes_over_stack(d00cm, {2,3}, 'color', 'y');

d00cm = stack_align(d00cm, {13,1}, {11,8});
show_electrodes_over_stack(d00cm, {13,1}, 'color', 'y');


%% Reslice
stack = {12,6};
d00cm.stacks{stack{:}} = load_slices(d00cm.stacks{stack{:}});
stack_reslicer(d00cm.stacks{stack{:}}, 'points_color', 'b', 'bgpoint_color', 'g', 'method', 'linear');
d00cm.stacks{stack{:}}.reslice_points = gui_points;
save_stacks(d00cm);

impath = '/confocal/2013-08-06/tile63x_17_18_resliced.tif';
imwrite(gui_resliced, fullfile(server_data_path(), piece, impath));


%% GUI mark EI soma and initial segment
phlmarks = fullfile(server_path(), piece, 'd00-01-02-r/data000', 'phlmarks');
load(phlmarks);

celltype = 1;
cellids = d00cm.cell_types{celltype}.cell_ids;
numcells = length(cellids);
for i = 1:numcells
    cellid = cellids(i);
    [~, marks1{i}] = mark_ei_somseg(d00cm, cellid);
end
backup = marks1;

celltype = 2;
cellids = d00cm.cell_types{celltype}.cell_ids;
numcells = length(cellids);
for i = 1:numcells
    cellid = cellids(i);
    [~, marks2{i}] = mark_ei_somseg(d00cm, cellid);
end
backup = marks2;

save(phlmarks, '-append', 'marks2');



%% Plot marks over anatomy!
ai = load_array_info(d00cm, 2);
stack = {13,1};

array_tform_inv = d00cm.stacks{stack{:}}.tforms_inv{1,1};
resliced = get_slice(d00cm.stacks, {stack{:}});
resliced(:,:,3) = resliced(:,:,2);
resliced(:,:,2) = resliced(:,:,1);
imshow(resliced); hold on

phlmarks = fullfile(server_path(), piece, 'd00-01-02-r/data000', 'phlmarks');
load(phlmarks);

for i = 1:length(marks1)
    mark = marks1{i};
    mark = tforminv(ai.T_array_to_vision_ei, mark);
    arrow = tformfwd(array_tform_inv, mark);

    h = imarrow(gca, arrow);
    c = get(h, 'Children');
    set(c, 'Color', 'w');
    set(c(4), 'Visible', 'off'); % Kill drop shadow
    
    mid = mean(arrow);
    h = text(mid(1), mid(2), sprintf('+P %d', d00cm.cell_types{1}.cell_ids(i)));
    set(h, 'Color', 'w');
    
    drawnow;
end

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
    h = text(mid(1), mid(2), sprintf('-P %d', d00cm.cell_types{2}.cell_ids(i)));
    set(h, 'Color', 'w');
    
    drawnow;
end



%% Save out nice EI example for GRI
cid = 1221;
frames = 1:30;

figure();
F = struct('cdata', [], 'colormap', []);
for f = frames
    cla;
    h = plot_ei(d00cm, cid, 'frame_number', f, 'pos_color', 'b', 'neg_color', 'b');
    plot_electrodes(d00cm.ei.position, 'scale', 0.1, 'pos_color', [0.5 0.5 0.5]);
    set(h(h>0), 'ZData', 100*ones(127,1))
    set(gca, 'Box', 'off', 'YTick', [], 'XTick', [], 'YColor', [1 1 1], 'XColor', [1 1 1])
    drawnow();
    F(end+1) = getframe();
end
close();
F = F(2:end);

vw = VideoWriter(sprintf('ei_%s_d00cm_%d', piece, cid));
vw.open();
vw.writeVideo(F);
vw.close();

% Projected frame
h = plot_ei(d00cm, cid, 'pos_color', 'b', 'neg_color', 'b');
plot_electrodes(d00cm.ei.position, 'scale', 0.1, 'pos_color', [0.5 0.5 0.5]);
set(h(h>0), 'ZData', 100*ones(127,1))
set(gca, 'Box', 'off', 'YTick', [], 'XTick', [], 'YColor', [1 1 1], 'XColor', [1 1 1])

