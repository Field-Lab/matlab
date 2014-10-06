piece = '2011-07-05-4';
d03 = load_data([piece '/data003']);
d03.piece.array_id = 1504;
d03.piece.rig = 'A';
d03 = restore_stacks(d03);
d03 = load_params(d03);
d03 = load_ei(d03, []);


%%
% open contourmatch_2011-07-05-4


%%
% d09 = load_data('2011-07-05-4/data009');
% d09.piece.array_id = 1504;
% d09.piece.rig = 'A';
% d09 = restore_stacks(d09);
% d09 = load_params(d09);
% d09 = load_ei(d09, []);


%% Should no longer need to be run
% datarun = init_stacks(datarun);
% datarun.stacks{2,1}  = build_stack('2011-07-05-4/spot/10xarrayb.tif');
% datarun.stacks{10,1} = build_stack('2011-07-05-4/confocal/2011-07-19/2011-07-05-4/2011-07-05-4_Series120_Merged_stack.tif:[1-14]');
% datarun.stacks{10,2} = build_stack('2011-07-05-4/confocal/2011-07-19/2011-07-05-4/2011-07-05-4_Series120_resliced1.tif');
% datarun.stacks{11,1} = build_stack('2011-07-05-4/confocal/2011-07-19/2011-07-05-4/2011-07-05-4_Series063_Merged_stack_RB.tif:[1-20]');
% datarun.stacks{11,2} = build_stack('2011-07-05-4/confocal/2011-07-19/2011-07-05-4/2011-07-05-4_Series063_Merged_stack.tif:[1-20]');
% datarun.stacks{11,3} = build_stack('2011-07-05-4/confocal/2011-07-19/2011-07-05-4/2011-07-05-4_Series063_RB_resliced1.tif');
% datarun.stacks{2,1}.basepath  = '/Users/vision/Desktop';
% datarun.stacks{10,1}.basepath  = '/Users/vision/Desktop';
% datarun.stacks{10,2}.basepath  = '/Users/vision/Desktop';
% datarun.stacks{11,1}.basepath = '/Users/vision/Desktop';
% datarun.stacks{11,2}.basepath = '/Users/vision/Desktop';
% datarun.stacks{11,3}.basepath = '/Users/vision/Desktop';
% 
% d03.stacks = fix_basepaths(d03.stacks);
% d03.stacks{11,10} = build_stack('2011-07-05-4/confocal/2011-07-19/2011-07-05-4/2011-07-05-4_Series048_062_Merge_stack.tif:[1-17]');
% d03.stacks{11,11} = build_stack('2011-07-05-4/confocal/2011-07-19/2011-07-05-4/2011-07-05-4_Series048_062_Merge_stack_resliced.tif');
% save_stacks(d03);

%imshow(get_slice(d03.stacks{11,8}))


%%
% load /snle/analysis/2011-07-05-4/Series063_reslice_points
% load /snle/analysis/2011-07-05-4/Series063_reslice_points2
% load /snle/analysis/2011-07-05-4/Series063_reslice_points3
% load /snle/analysis/2011-07-05-4/Series120_reslice_points
% load /snle/data/2011-07-05-4/confocal/2011-07-19/2011-07-05-4/gui_points.mat

% Move this into stacks!
% load(fullfile(server_path, piece, 'Series063_reslice_points4'));
% d03.stacks{11,1}.reslice_points = Series063_reslice_points4;
% d03.stacks{11,2}.reslice_points = Series063_reslice_points4;
% save_stacks(d03)

%%
% resliced3 = get_slice(d03.stacks{11,3});
% imshow(resliced3(:,:,1));
% plot_parasol_contours(d03);


%% Reslice
d03.stacks{11,10} = load_slices(d03.stacks{11,10});
stack_reslicer(d03.stacks{11,10}, 'points_color', 'b', 'bgpoint_color', 'g', 'method', 'linear', 'channels', 1);

d03.stacks{11,10}.reslice_points = gui_points;
save_stacks(d03)

[p,n,e] = fileparts(d03.stacks{11,10}.paths{1});
imwrite(gui_resliced, fullfile(d03.stacks{11,10}.basepath(), p, sprintf('%s_resliced%s', n, e)));
clear p n e;


%% Mark anatomical cell centers
stack_point_picker(d03.stacks{11,11}, 'points_color', 'w', 'bgpoint_color', 'b', 'channels', 1, 'points', {anatmarks2});
anatmarks2 = gui_points{1};
save(fullfile(server_path(), piece, 'anatmarks'), 'anatmarks1', 'anatmarks2')

anatmarks = anatmarks1;
for i = 1:length(anatmarks)
    plot(anatmarks(i,1), anatmarks(i,2), 'rs');
end; clear i anatmarks


%% HACKETY fix for horrible point duplication issue
% new_gui_points = cell(1,length(gui_points));
% 
% for i = 1:length(gui_points)
%     points = gui_points{i};
%     if isempty(points), continue; end
%     
%     [u,n,m] = unique(points(:,1));
%     new_gui_points{i} = points(n,:);
% end
% 
% clear u n m i points


%%
% stack_align(d03, {11,9}, {10,2});

%% Plot ei over image

cellid = d03.cell_types{1}.cell_ids(2);

plot_ei_scroll(d03, cellid, 'coordinates', 'stack', 'stack', {11,9}, 'max_scale', 0.75, 'scale', 1.5, 'pos_color', 'g');
eiax = gca;
set(eiax, 'Color', 'none');

stackax = axes;
imshow(resliced(:,:,1));
stackax = gca;

linkaxes([eiax stackax]);
axes(eiax);

set(gcf, 'Units', 'normalized', 'Position', [0.5 0 0.5 1]);
autozoom_ei(d03, cellid, 'stack', {11,3});
axis fill;



%% Plot current source densities
celltype = 1;

ei_modification = @(ei,datarun,cellid)(ei2csd(ei, d03.ei.neighbor_struct));
for i = 1:length(d03.cell_types{celltype}.cell_ids)
    cellid = d03.cell_types{celltype}.cell_ids(i);
    plot_ei_scroll(d03, cellid, 'modify_ei', ei_modification, 'max_scale', 0.75, 'pos_color', 'g', 'figure', i);
    set(gca, 'Color', 'none');
end

% Nice ON Parasols: 2 5 34
% Nice OFF Parasols: 5 9 26
% Nice ON Midgets: 6 11 16 36 51 66 71


%% Try to model fit EI or CSD
% First attempt: just two round 2D gaussians with gaussian timecourse


%% GUI locate soma and initial segment
d03 = load_ei(d03, {4});
d03 = calc_csd(d03, 'cellspec', {4});

% ON Parasols
celltype = 1;
cellids = d03.cell_types{celltype}.cell_ids;
numcells = length(cellids);
load(fullfile(server_path, piece, 'data003/marks1'));
for i = 22
    cellid = cellids(i);
    [~, marks1{i}] = mark_ei_somseg(d03, cellid, 'marks', marks1{i});
end
backup = marks1;
save /snle/analysis/2011-07-05-4/data003/marks1 marks1


% OFF Parasols
celltype = 2;
cellids = d03.cell_types{celltype}.cell_ids;
numcells = length(cellids);
load(fullfile(server_path, piece, 'data003/marks2'));

for i = numcells:numcells
    cellid = cellids(i);
    [~, marks2{i}] = mark_ei_somseg(d03, cellid, 'marks', marks2{i});
end
backup = marks2;

% ON Midgets
celltype = 3;
cellids = d03.cell_types{celltype}.cell_ids;
numcells = length(cellids);
for i = 1:numcells
    cellid = cellids(i);
    [~, marks3{i}] = mark_ei_somseg(d03, cellid);
end
backup = marks3;

% OFF Midgets
celltype = 4;
cellids = d03.cell_types{celltype}.cell_ids;
numcells = length(cellids);
for i = length(marks4):numcells
    cellid = cellids(i);
    if length(marks4) >= i
        [~, marks4{i}] = mark_ei_somseg(d03, cellid, 'marks', marks4{i});
    else
        [~, marks4{i}] = mark_ei_somseg(d03, cellid);
    end
end
backup = marks4;

% Amacrines
celltype = 6;
cellids = d03.cell_types{celltype}.cell_ids;
numcells = length(cellids);
for i = 1:(numcells-1)
    cellid = cellids(i);
    [~, marks6{i}] = mark_ei_somseg(d03, cellid);%, 'marks', marks6{i});
end
backup = marks6;

celltype = 16;
cellids = d03.cell_types{celltype}.cell_ids;
numcells = length(cellids);
for i = 1:(numcells-1)
    cellid = cellids(i);
    [~, marks16{i}] = mark_ei_somseg(d03, cellid);%, 'marks', marks6{i});
end


%% PHL revised marks
load(fullfile(server_path, piece, 'data003/marksPHL'));
save(fullfile(server_path, piece, 'data003/marksPHL'), 'marks1', 'marks2', 'marks3', 'marks4', 'marks6', 'marks16');


%% Plot marks over anatomy!
datarun = d03;
ai = load_array_info(datarun, 2);
stack = {11,11};

datarun.stacks{stack{:}}.basepath = server_data_path();
array_tform_inv = datarun.stacks{stack{:}}.tforms_inv{1,1};
resliced = get_slice(datarun.stacks, {stack{:}});
resliced = resliced(:,:,1);
imshow(resliced); hold on

% % Just soma mode
% for i = 1:length(marks1)
%     mark = marks1{i};
%     mark = tforminv(ai.T_array_to_vision_ei, mark);
%     soma = tformfwd(array_tform_inv, mark(1,:));
% 
%     h = plot(soma(1), soma(2), 'o');
%     set(h, 'Color', 'b', 'LineWidth', 2);
%     
% %     h = text(soma(1)+5, soma(2)-5, sprintf('+P %d', datarun.cell_types{1}.cell_ids(i)), 'FontSize', 36);
% %     set(h, 'Color', 'b');
% end
% drawnow;
% 
% for i = 1:length(marks2)
%     mark = marks2{i};
%     mark = tforminv(ai.T_array_to_vision_ei, mark);
%     soma = tformfwd(array_tform_inv, mark(1,:));
% 
%     h = plot(soma(1), soma(2), 'o');
%     set(h, 'Color', 'g', 'LineWidth', 2);
%     
% %     h = text(soma(1)+5, soma(2)-5, sprintf('-P %d', datarun.cell_types{2}.cell_ids(i)), 'FontSize', 36);
% %     set(h, 'Color', 'g');
% end
% drawnow;

% Arrow mode
for i = 1:length(marks1)
    mark = marks1{i};
    mark = tforminv(ai.T_array_to_vision_ei, mark);
    arrow = tformfwd(array_tform_inv, mark);

    h = imarrow(gca, arrow);
    c = get(h, 'Children');
    set(c, 'Color', 'y');
    set(c(4), 'Visible', 'off'); % Kill drop shadow
    
%     mid = mean(arrow);
%     h = text(mid(1), mid(2), sprintf('+P %d', datarun.cell_types{1}.cell_ids(i)));
%     set(h, 'Color', 'b');
    
    drawnow;
end
for i = 1:length(marks2)
    mark = marks2{i};
    mark = tforminv(ai.T_array_to_vision_ei, mark);
    arrow = tformfwd(array_tform_inv, mark);

    h = imarrow(gca, arrow);
    c = get(h, 'Children');
    set(c, 'Color', 'b');
    set(c(4), 'Visible', 'off'); % Kill drop shadow
    
%     mid = mean(arrow);
%     h = text(mid(1), mid(2), sprintf('+P %d', datarun.cell_types{1}.cell_ids(i)));
%     set(h, 'Color', 'b');
    
    drawnow;
end


%% Permutation test: average distance from anatmark to nearest mark
% anatmarks windowed, marks randomly shifted, rotated, flipped
load(fullfile(server_path, piece, 'data003/permtest'))
load(fullfile(server_path, piece, 'data003/marksPHL'));
load(fullfile(server_path(), piece, 'anatmarks'))

anatmarks = anatmarks2;
eimarks = marks4;
saveindex = {2,4}; % Should match anatmarks, eimarks
datarun = d03;
ai = load_array_info(datarun, 2);

% FIXME: Should use 11,11 but apparently that alignment still needs to be
% synced from my home computer :(
stack = {11,11};
% stack = {11,9};

array_tform_inv = datarun.stacks{stack{:}}.tforms_inv{1,1};

somas = cell2mat(cellfun(@(M) M(1,:), eimarks, 'UniformOutput', false)');
somas = tforminv(ai.T_array_to_vision_ei, somas);
somas = tformfwd(array_tform_inv, somas);

window = [650 2650 1300 3300];
anatmarks = anatmarks(anatmarks(:,1) > window(1),:);
anatmarks = anatmarks(anatmarks(:,1) < window(2),:);
anatmarks = anatmarks(anatmarks(:,2) > window(3),:);
anatmarks = anatmarks(anatmarks(:,2) < window(4),:);

d = ipdm(anatmarks, somas, 'Subset', 'NearestNeighbor');
d = d(d ~= 0);
m = mean(d);

somhull = convhull(somas(:,1), somas(:,2));
window_center = [window(1)+window(2) window(3)+window(4)] ./ 2;
window_range = [window(2) - window(1), window(4) - window(3)];

t_soma_collect = {somas};
d_collect = full(d);

for i = 1:10000
    % Transform soma locations
    t_somas = somas;
    
    % Rotate
    rdeg = randi([1 360]);
    rmat = [cosd(rdeg) sind(rdeg); -sind(rdeg) cosd(rdeg)];
    t_somas = t_somas*rmat';
    
    % Move to origin
    soma_center = centroid(t_somas(somhull,1), t_somas(somhull,2));
    t_somas = t_somas - repmat(soma_center, size(t_somas,1), 1);
    
    % Flip?
    if randi([0 1]), t_somas(:,1) = -t_somas(:,1); end
    
    % Shift
    shiftrange = (max(t_somas) - min(t_somas) - window_range) ./ 2.5; % heuristic scale factor
    xshift = randi(round([-shiftrange(1) shiftrange(1)]));
    yshift = randi(round([-shiftrange(2) shiftrange(2)]));
    t_somas = t_somas + repmat([xshift yshift], size(t_somas,1), 1);
    
    % Move to window
    t_somas = t_somas + repmat(window_center, size(t_somas,1), 1);
    
    
    % Test whether window is still inside convhull after shifting somas
    if ~all(inpolygon(window([1; 1; 2; 2]), window([3; 4; 3; 4]), t_somas(somhull,1), t_somas(somhull,2)))
        continue;
    end
    
    % Test that points have moved sufficiently
    d = sqrt(sum((somas - t_somas).^2, 2));
    if min(d) <= 250 % Roughly based on point spacing
        continue;
    end
    
    d = ipdm(anatmarks, t_somas, 'Subset', 'NearestNeighbor');
    t_soma_collect{end+1} = t_somas;
    d_collect(:,end+1) = d(d > 0);

%     % Visualize
%     cla
%     hold on
%     plot(somas(:,1), somas(:,2));
%     plot(t_somas(:,1), t_somas(:,2), 'g');
%     axis equal
%     r = rectangle('Position', [window(1), window(3), window(2)-window(1), window(4)-window(3)], 'EdgeColor', 'b');
    
end; clear i xshift yshift shiftrange t_somas soma_center rmat rdeg
hist(mean(d_collect), 100)
drawnow

permd{saveindex{:}} = d_collect;
permsom{saveindex{:}} = t_soma_collect;
trued(saveindex{:}) = m;
truesom{saveindex{:}} = somas;
save(fullfile(server_path, piece, 'data003/permtest'), 'permd', 'permsom', 'trued', 'truesom')
clear d_collect t_soma_collect somas d m

%% Plot permutation test results
% load(fullfile(server_path, piece, 'data003/permtest'))

figure();
set(gcf, 'Color', 'w', 'Units', 'normalized', 'Position', [0.75 0 0.25 1])

%%
clf
ylims = {[0 7000] [0 7000] [0 6000] [0 12000]};

dx = 5;
histx = 0:dx:450;
for r = 1:size(permd,1)
    for c = 1:size(permd,2)
        sanesubplot(size(permd,2), size(permd,1), {c, r});
        histy = hist(mean(permd{r,c}), histx);
        b = bar(histx, histy, 'hist');
        set(b, 'FaceColor', 'none');
        set(gca, 'YLim', ylims{c});
        set(gca, 'XLim', [0 450]);
        set(gca, 'Color', 'none', 'Box', 'off', 'YTick', [], 'YColor', 'w');
        
        % These annotations are annoying.  Only works if you get figure
        % window to desired size before calculating annotation location...
        [figx figy] = dsxy2figxy([trued(r,c), trued(r,c)], [0, 0]);
        annotation('arrow', figx, figy - [0.05 0.01]);
    end
end

% X Labels (drawn manually since we set YColor to white to hide the axis line
sanesubplot(4,2, {1 1});
text(100, 8000, 'Large, brightly HCN1 immunoreactive, polygonal somas');
sanesubplot(4,2, {1 2});
text(100, 8000, 'Round, dimly HCN1 immunoreactive somas');

% Y Labels
sanesubplot(4,2, {1 1});
text(-20, 200, 'ON Parasol electrical images', 'rotation', 90);
sanesubplot(4,2, {2 1});
text(-20, 200, 'OFF Parasol electrical images', 'rotation', 90);
sanesubplot(4,2, {3 1});
text(-20, 200, 'ON Midget electrical images', 'rotation', 90);
sanesubplot(4,2, {4 1});
text(-20, 200, 'OFF Midget electrical images', 'rotation', 90);


%% Compare EI and CSD
celltype = 1;
xlims = [10 25];

for i = 1:length(d03.cell_types{celltype}.cell_ids)
    cellid = d03.cell_types{celltype}.cell_ids(i);
    ei = get_ei(d03, cellid);
    csd = ei2csd(ei, d03.ei.neighbor_struct);

    figure
    sanesubplot(1,2,1);
    plot(ei');
    set(gca, 'XLim', xlims);
    grid on;
    
    sanesubplot(1,2,2);
    plot(csd');
    set(gca, 'XLim', xlims);
    grid on;
end



%% Plot current source density over image

cellid = d03.cell_types{1}.cell_ids(2);

ei_modification = @(ei,datarun,cellid)(ei2csd(ei, d03.ei.neighbor_struct));
plot_ei_scroll(d03, cellid, 'coordinates', 'stack', 'stack', {11,9}, 'modify_ei', ei_modification, 'max_scale', 0.75, 'scale', 3, 'pos_color', 'g');
eiax = gca;
set(eiax, 'Color', 'none');

stackax = axes;
imshow(resliced(:,:,1));
stackax = gca;

linkaxes([eiax stackax]);
axes(eiax);

set(gcf, 'Units', 'normalized', 'Position', [0.56 0 0.45 1]);
%autozoom_ei(d03, cellid, 'stack', {11,3}, 'padding_factor', 5, 'keep_aspect_ratio', false);
%axis fill;
