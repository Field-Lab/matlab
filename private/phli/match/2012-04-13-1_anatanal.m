%%
piece = '2012-04-13-1';
rig = 'A';
loadopts = struct('load_params', true, 'load_neurons', true, 'load_ei', true);
keep_vars = {'piece'; 'loadopts'; 'keep_vars'};


%%
d00 = load_data([piece '/data000'], loadopts);
d00 = restore_stacks(d00);

datarun = d00;

keep_vars = [keep_vars; 'd00'; 'datarun'];

% d03.piece.array_id = 1504;
% d03.piece.rig = 'A';
% d03 = restore_stacks(d03);

%% EI Marking
celltype = 6;
cellids = d00.cell_types{celltype}.cell_ids;
numcells = length(cellids);
load(fullfile(server_path, piece, sprintf('data000/marks%d', celltype)));

for i = 1:numcells
    cellid = cellids(i);
    [~, marks6{i}] = mark_ei_somseg(d00, cellid);
end
backup = marks6;
save(fullfile(server_path, piece, sprintf('data000/marks%d', celltype)), sprintf('marks%d', celltype));

%% Plot marks over anatomy!
ai = load_array_info(d00, 2);
stackind = {11,2};

d00.stacks{stackind{:}}.basepath = server_data_path();
array_tform_inv = d00.stacks{stackind{:}}.tforms_inv{1,1};
resliced = get_slice(d00.stacks, {stackind{:}});
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

load(fullfile(server_path, piece, 'data000/marks4'));
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
    h = text(mid(1), mid(2), sprintf('-M %d', d00.cell_types{4}.cell_ids(i)));
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

%% Setting up stacks; only run once then use restore_stacks
% datarun = init_stacks(datarun);
% datarun.stacks{2,1}  = build_stack([piece '/spot/6.5x/stitched.tif']);
% datarun.stacks{2,2}  = build_stack([piece '/spot/15x/stitched.tif']);
%
% datarun = array_align(datarun, {2,1}, 'rig', rig);
% datarun = array_align(datarun, {2,2}, 'rig', rig);
% 
% datarun.stacks{10,1} = build_stack([piece '/confocal/2012-04-13-1_Brn3_NHP/Brn3_NHP_Series012_Merged_stack.tif:[1-8]']);
% datarun.stacks{10,2} = build_stack([piece '/confocal/2012-04-13-1_Brn3_NHP/Brn3_NHP_Series047_Merged_stack.tif:[1-16]']);
% datarun.stacks{11,1} = build_stack([piece '/confocal/2012-04-13-1_Brn3_NHP/Brn3_NHP_Series012_Merged_stack_resliced.tif']);
% datarun.stacks{11,2} = build_stack([piece '/confocal/2012-04-13-1_Brn3_NHP/Brn3_NHP_Series047_Merged_stack_resliced.tif']);
% d00.stacks = fix_basepaths(d00.stacks);
% d00.stacks{10,3} = build_stack(fullfile(piece, '/confocal/2012-04-13-1_Brn3_DAPI/2012-04-13-1_Brn3_DAPI_stitch.tif:[1-6]'));
%
% save_stacks(d00)
%  
% figure; imshow(get_slice(datarun.stacks{11,2}))


%% Reslice stack
reslicestack = {10,3};

% Pick points
d00.stacks{reslicestack{:}} = load_slices(d00.stacks{reslicestack{:}});
stack_reslicer(d00.stacks{reslicestack{:}}, 'points_color', 'b', 'bgpoint_color', 'r', 'method', 'linear');
d00.stacks{reslicestack{:}}.reslice_points = gui_points;

%% Save reslice points
% Generally the stack_reslicer GUI saves points to the base workspace with
% name gui_points.  Kind of kludgy but good enough for now...
% datarun.stacks{reslicestack{:}}.reslice_points = gui_points;
% save_stacks(datarun);


%% Save resliced image

% resliced = getappdata(gcf, 'resliced');
% [p,f] = fileparts(fullfile(datarun.stacks{reslicestack{:}}.basepath, datarun.stacks{reslicestack{:}}.paths{1}));
% imwrite(resliced, sprintf('%s/%s%s.tif', p, f, '_resliced'));

% Now loop back up and add to stacks


%% stack_align
% datarun = stack_align(datarun, {11,1}, {2,2});
% datarun = stack_align(datarun, {11,2}, {11,1});


%% Build tforms
% stacks = build_tform_routes(stacks, varargin)
% datarun.stacks = build_tform_from_routes(datarun.stacks, {11,1}, {1,1});
% datarun.stacks = build_tform_from_routes(datarun.stacks, {11,2}, {1,1});


%% Sanity check
figure;show_electrodes_over_stack(datarun, {11,1}, 'color', 'g', 'include_outside', true);


%% Plot just max electrode location
image = get_slice(datarun.stacks{11,2});
figure; imshow(image(:,:,1));
hold on
array_tform_inv = datarun.stacks{11,2}.tforms_inv{1,1};

cell_type = 1;
for i = 1:length(datarun.cell_types{cell_type}.cell_ids)
    id = datarun.cell_types{cell_type}.cell_ids(i);
    maxelectrode = datarun.ei.java_ei.getMaxElectrode(id);
    [x,y] = tformfwd(array_tform_inv, datarun.ei.position(maxelectrode,:));
    h = plot(x,y,'bo');
    set(h, 'LineWidth', 2)
end


%% Prep for plot physiology over image
datarun = calc_csd(datarun);




%% Below just copied over from 2011-07-05-4; need to convert everything over.


%% Plot ei over image

% cellid = d03.cell_types{1}.cell_ids(2);
% 
% plot_ei_scroll(d03, cellid, 'coordinates', 'stack', 'stack', {11,9}, 'max_scale', 0.75, 'scale', 1.5, 'pos_color', 'g');
% eiax = gca;
% set(eiax, 'Color', 'none');
% 
% stackax = axes;
% imshow(resliced(:,:,1));
% stackax = gca;
% 
% linkaxes([eiax stackax]);
% axes(eiax);
% 
% set(gcf, 'Units', 'normalized', 'Position', [0.5 0 0.5 1]);
% autozoom_ei(d03, cellid, 'stack', {11,3});
% axis fill;
% 


%% Plot current source densities
celltype = 1;
 
ei_modification = @(ei,datarun,cellid)(ei2csd(ei, datarun.ei.neighbor_struct));
for i = 1:length(datarun.cell_types{celltype}.cell_ids)
    cellid = datarun.cell_types{celltype}.cell_ids(i);
    plot_ei_scroll(datarun, cellid, 'modify_ei', ei_modification, 'max_scale', 0.75, 'pos_color', 'g', 'figure', i);
    set(gca, 'Color', 'none');
end


%% Compare EI and CSD
% celltype = 1;
% xlims = [10 25];
% 
% for i = 1:length(d03.cell_types{celltype}.cell_ids)
%     cellid = d03.cell_types{celltype}.cell_ids(i);
%     ei = get_ei(d03, cellid);
%     csd = ei2csd(ei, d03.ei.neighbor_struct);
% 
%     figure
%     sanesubplot(1,2,1);
%     plot(ei');
%     set(gca, 'XLim', xlims);
%     grid on;
%     
%     sanesubplot(1,2,2);
%     plot(csd');
%     set(gca, 'XLim', xlims);
%     grid on;
% end



%% Plot current source density over image

cellid = datarun.cell_types{1}.cell_ids(19);
 
ei_modification = [];%@(ei,datarun,cellid)(ei2csd(ei, datarun.ei.neighbor_struct));
plot_ei_scroll(datarun, cellid, 'coordinates', 'stack', 'stack', {11,2}, 'modify_ei', ei_modification, 'max_scale', 0.75, 'scale', 3, 'pos_color', 'g');
eiax = gca;
set(eiax, 'Color', 'none');
 
stackax = axes;
imshow(image(:,:,1));
stackax = gca;
 
linkaxes([eiax stackax]);
axes(eiax);
 
set(gcf, 'Units', 'normalized', 'Position', [0.56 0 0.45 1]);
% autozoom_ei(datarun, cellid, 'stack', {11,3}, 'padding_factor', 5, 'keep_aspect_ratio', false);
% axis fill;
