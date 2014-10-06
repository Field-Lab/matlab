%%
%load contourmatch_2011-07-05-4.mat
load 2011-07-05-4_maxcontours.mat

%%
piece = '2011-07-05-4';
d03 = load_data('2011-07-05-4/data003');
d03.piece.array_id = 1504;
d03.piece.rig = 'A';
d03 = restore_stacks(d03);
d03 = load_params(d03);
d03 = load_ei(d03, []);


%% Use reslice points from RB to also reslice RGB
% load /snle/data/2011-07-05-4/confocal/2011-07-19/2011-07-05-4/Series063_reslice_points2.mat
% d03.stacks{11,2} = load_slices(d03.stacks{11,2});
% stack_reslicer(d03.stacks{11,2}, 'points', Series063_reslice_points2, 'method', 'linear');
% % Add points if you want, pick "reslice stack" from GUI menu
% resliced = getappdata(gcf, 'resliced');

% Now just stored in stacks{11,9}


%% {11,3}
array_tform_inv = d03.stacks{11,3}.tforms_inv{1,1};
resliced3 = get_slice(d03.stacks, {11,3});
imshow(resliced3(:,:,1)); hold on


%% {11,5}
array_tform_inv = d03.stacks{11,5}.tforms_inv{1,1};
resliced = get_slice(d03.stacks, {11,5});
imshow(resliced); hold on


%% {11,7}
array_tform_inv = d03.stacks{11,7}.tforms_inv{1,1};
resliced = get_slice(d03.stacks, {11,7});
imshow(resliced(:,:,1)); hold on


%% {11,9}
array_tform_inv = d03.stacks{11,9}.tforms_inv{1,1};
resliced = get_slice(d03.stacks, {11,9});
imshow(resliced(:,:,1)); hold on


%% DPA EI soma and initial segment marks
load(fullfile(server_path, piece, 'data003/marks1DA'));
load(fullfile(server_path, piece, 'data003/marks2DA'));
ai = load_array_info(d03, 2);

for i = 1:length(marks1)
    mark = marks1{i};
    mark = tforminv(ai.T_array_to_vision_ei, mark);
    arrow = tformfwd(array_tform_inv, mark);
    h = imarrow(gca, arrow);
    c = get(h, 'Children');
    set(c, 'Color', [1 0 1]);
    set(c(4), 'Visible', 'off'); % Turn off drop shadow
    drawnow;
end

for i = 1:length(marks2)
    mark = marks2{i};
    mark = tforminv(ai.T_array_to_vision_ei, mark);
    arrow = tformfwd(array_tform_inv, mark);
    h = imarrow(gca, arrow);
    c = get(h, 'Children');
    set(c, 'Color', 'g');
    set(c(4), 'Visible', 'off'); % Turn off drop shadow
    drawnow;
end


%% PHL revised marks
load(fullfile(server_path, piece, 'data003/marksPHL'));
ai = load_array_info(d03, 2);

for i = 1:length(marks1)
    mark = marks1{i};
    mark = tforminv(ai.T_array_to_vision_ei, mark);
    arrow = tformfwd(array_tform_inv, mark);
    h = imarrow(gca, arrow);
    c = get(h, 'Children');
    set(c, 'Color', [0 0 1]);
    set(c(4), 'Visible', 'off'); % Turn off drop shadow
    drawnow;
end

for i = 1:length(marks2)
    mark = marks2{i};
    mark = tforminv(ai.T_array_to_vision_ei, mark);
    arrow = tformfwd(array_tform_inv, mark);
    h = imarrow(gca, arrow);
    c = get(h, 'Children');
    set(c, 'Color', 'g');
    set(c(4), 'Visible', 'off'); % Turn off drop shadow
    drawnow;
end


%% Check amacrines
% Pretty hard to tell, but does seem like there is some correspondence
% between amacrine physiological marks and "funny nuclei, thin soma" cells.

array_tform_inv = d03.stacks{11,9}.tforms_inv{1,1};
resliced = get_slice(d03.stacks, {11,9});
imshow(resliced); hold on

load(fullfile(server_path, piece, 'data003/marksPHL'));
ai = load_array_info(d03, 2);
for i = 1:length(marks1)
    mark = marks1{i};
    mark = tforminv(ai.T_array_to_vision_ei, mark);
    arrow = tformfwd(array_tform_inv, mark);
    
    h = imarrow(gca, arrow);
    c = get(h, 'Children');
    set(c, 'Color', [1 1 1]);
    set(c(3), 'LineStyle', '--');
    set(c(4), 'Visible', 'off'); % Turn off drop shadow

    mid = mean(arrow);
    h = text(mid(1), mid(2), sprintf('P+%d', d03.cell_types{1}.cell_ids(i)));
    set(h, 'Color', 'w');

    drawnow;
end
for i = 1:length(marks2)
    mark = marks2{i};
    mark = tforminv(ai.T_array_to_vision_ei, mark);
    arrow = tformfwd(array_tform_inv, mark);

    h = imarrow(gca, arrow);
    c = get(h, 'Children');
    set(c, 'Color', [1 1 1]);
    set(c(3), 'LineStyle', '-.');
    set(c(4), 'Visible', 'off'); % Turn off drop shadow

    mid = mean(arrow);
    h = text(mid(1), mid(2), sprintf('P-%d', d03.cell_types{2}.cell_ids(i)));
    set(h, 'Color', 'w');

    drawnow;
end
for i = 1:length(marks6)
    mark = marks6{i};
    mark = tforminv(ai.T_array_to_vision_ei, mark);
    arrow = tformfwd(array_tform_inv, mark);

    h = imarrow(gca, arrow);
    c = get(h, 'Children');
    set(c, 'Color', [1 1 1]);
    set(c(3), 'LineStyle', '-.');
    set(c(4), 'Visible', 'off'); % Turn off drop shadow

    mid = mean(arrow);
    h = text(mid(1), mid(2), sprintf('A-%d', d03.cell_types{6}.cell_ids(i)));
    set(h, 'Color', 'w');

    drawnow;
end
for i = 1:length(marks16)
    mark = marks16{i};
    mark = tforminv(ai.T_array_to_vision_ei, mark);
    arrow = tformfwd(array_tform_inv, mark);

    h = imarrow(gca, arrow);
    c = get(h, 'Children');
    set(c, 'Color', [1 1 1]);
    set(c(4), 'Visible', 'off'); % Turn off drop shadow

    mid = mean(arrow);
    h = text(mid(1), mid(2), sprintf('A+%d', d03.cell_types{16}.cell_ids(i)));
    set(h, 'Color', 'w');

    drawnow;
end


%% Cell nums at max electrode locations
for i = 1:length(d03.cell_types{1}.cell_ids)
    id = d03.cell_types{1}.cell_ids(i);
    maxelectrode = d03.ei.java_ei.getMaxElectrode(id);
    [x,y] = tformfwd(array_tform_inv, d03.ei.position(maxelectrode,:));
    h = text(x,y,num2str(i), 'Color', 'b');
end

for i = 1:length(d03.cell_types{2}.cell_ids)
    id = d03.cell_types{2}.cell_ids(i);
    maxelectrode = d03.ei.java_ei.getMaxElectrode(id);
    [x,y] = tformfwd(array_tform_inv, d03.ei.position(maxelectrode,:));
    h = text(x,y,num2str(i), 'Color', 'g');
end


%% Max contours
for i = 1:length(d03.cell_types{1}.cell_ids)
    cell_num = get_cell_indices(d03, d03.cell_types{1}.cell_ids(i));
    contour = maxcontours{cell_num};
    for j = length(contour)%:-10:1
        path = contour(j).path;
        path = tformfwd(array_tform_inv, path')';
        h = plot(path(1,:), path(2,:), 'b');
        set(h, 'LineWidth', 2);
        drawnow;
    end
end


for i = 1:length(d03.cell_types{2}.cell_ids)
    cell_num = get_cell_indices(d03, d03.cell_types{2}.cell_ids(i));
    contour = maxcontours{cell_num};
    for j = length(contour)%:-5:length(contour)-25
        path = contour(j).path;
        path = tformfwd(array_tform_inv, path')';
        h = plot(path(1,:), path(2,:), 'g');
        set(h, 'LineWidth', 1.5);
        drawnow;
    end
end


for i = 1:length(d03.cell_types{3}.cell_ids)
    cell_num = get_cell_indices(d03, d03.cell_types{3}.cell_ids(i));
    contour = maxcontours{cell_num};
    for j = length(contour) %[length(contour) length(contour)-3 length(contour)-6]
        path = contour(j).path;
        path = tformfwd(array_tform_inv, path')';
        h = plot(path(1,:), path(2,:), 'm-.');
        set(h, 'LineWidth', 1.5);
        drawnow;
    end
end


for i = 1:length(d03.cell_types{4}.cell_ids)
    cell_num = get_cell_indices(d03, d03.cell_types{4}.cell_ids(i));
    contour = maxcontours{cell_num};
    for j = length(contour) %[length(contour) length(contour)-3 length(contour)-6]
        path = contour(j).path;
        path = tformfwd(array_tform_inv, path')';
        h = plot(path(1,:), path(2,:), 'c-.');
        set(h, 'LineWidth', 1.5);
        drawnow;
    end
end


for i = 1:length(d03.cell_types{6}.cell_ids)
    cell_num = get_cell_indices(d03, d03.cell_types{6}.cell_ids(i));
    contour = maxcontours{cell_num};
    for j = length(contour) %[length(contour) length(contour)-3 length(contour)-6]
        path = contour(j).path;
        path = tformfwd(array_tform_inv, path')';
        h = plot(path(1,:), path(2,:), 'k-.');
        set(h, 'LineWidth', 1.5);
        drawnow;
    end
end

clear h i j


%% Compare to just max electrode location
for i = 1:length(d03.cell_types{1}.cell_ids)
    id = d03.cell_types{1}.cell_ids(i);
    maxelectrode = d03.ei.java_ei.getMaxElectrode(id);
    [x,y] = tformfwd(array_tform_inv, d03.ei.position(maxelectrode,:));
    h = plot(x,y,'bo');
    set(h, 'LineWidth', 2)
end

for i = 1:length(d03.cell_types{2}.cell_ids)
    id = d03.cell_types{2}.cell_ids(i);
    maxelectrode = d03.ei.java_ei.getMaxElectrode(id);
    [x,y] = tformfwd(array_tform_inv, d03.ei.position(maxelectrode,:));
    h = plot(x,y,'go');
    set(h, 'LineWidth', 2)
end


%%
d03 = load_ei(d03, []);
ei = get_ei_max_frame(get_ei(d03, d03.cell_types{1}.cell_ids(4)));
[cx,cy,coeffs] = compute_3dbspline_coeffs(d03.ei.position(:,1), d03.ei.position(:,2), abs(ei));
min(d03.ei.position(:,1))
min(d03.ei.position(:,2))
xi = -360:360;
yi = -390:390;
[XI,YI] = meshgrid(xi,yi);
ZI = interp_3dbspline(cx,cy,coeffs,XI,YI,2);


%% Plot surface
surf(XI, YI, ZI, 'EdgeColor', 'none');


%% Plot contours
numcontours = 50;
C = parse_contourc(xi, yi, ZI, numcontours);

hold on;
for i = 1:length(C)
    plot(C(i).path(1,:), C(i).path(2,:));
    drawnow;
end


%% Futzing...
contour = contours{1};
tic
for j = 1:length(contour)
    path = contour(j).path;
    h = plot3(path(1,:), path(2,:), repmat(contour(j).elevation, size(path(1,:))));
end
toc