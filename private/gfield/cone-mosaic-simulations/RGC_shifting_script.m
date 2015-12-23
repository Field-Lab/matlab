% this script shuffles the RGCs in a data set to the cone mosaic in the
% same dataset
clear all

cell_type =4;
radius_scaler = 2;
data_fruit = 'plantain';
jitter_scale = 4;
RGC_min_convergence = 7;
RGC_max_convergence = 27;

datarun = import_single_cone_data([], data_fruit);

% exclude S and U cones from type matrix, location matrix, connectivity
% matrix and roi matrix
S_cone_indices = find(datarun.cone_types == 'S');
U_cone_indices = find(datarun.cone_types == 'U');
exclude_cone_indices = [S_cone_indices; U_cone_indices];
keep_cone_indices = setdiff([1:length(datarun.cone_ids)], exclude_cone_indices);
datarun.cone_centers = datarun.cone_centers(keep_cone_indices,:);
datarun.cone_ids = datarun.cone_ids(keep_cone_indices);
datarun.cone_roi = datarun.cone_roi(keep_cone_indices);
datarun.cone_types = datarun.cone_types(keep_cone_indices);
datarun.cone_weights = datarun.cone_weights(keep_cone_indices,:);
datarun

RGC_indices = find(datarun.cell_types == cell_type);
num_RGCs = length(RGC_indices);
num_cones = length(datarun.cone_centers);


% precompute the new RGC locations
rand('twister', 1111);
location_jitter = jitter_scale * rand(num_RGCs,2);
roi_cones = find(datarun.cone_roi == 1);
num_roi_cones = length(roi_cones);
rand_cone_indices = randperm(num_roi_cones);
kosher_locations = datarun.cone_centers(roi_cones(rand_cone_indices(1:num_RGCs)),:) + location_jitter;

new_connectivity_matrix = zeros(num_cones, num_RGCs);
new_RGC_locations = zeros(num_RGCs,2);
excluded_cell_counter = 0;
for RGC = 1:length(RGC_indices)
    
    % get number of center cones and associated weights
    temp_RGC_center = datarun.rgc_COMs(RGC_indices(RGC),:);
    temp_RGC_CA = datarun.rgc_fit_info(RGC_indices(RGC),3) * radius_scaler;
    temp_distances = ipdm(temp_RGC_center, datarun.cone_centers);
    close_cone_indices = find(temp_distances <= temp_RGC_CA);
    positive_cone_indices = find(datarun.cone_weights(:, RGC_indices(RGC)) > 0);
    temp_cone_indices = intersect(close_cone_indices, positive_cone_indices);
    temp_num_center_cones = length(temp_cone_indices);
    
    clear temp_cone_weights
    if temp_num_center_cones < 7 || temp_num_center_cones > 27
        excluded_cell_counter = excluded_cell_counter + 1;
        temp_cone_weights(1:temp_num_center_cones) = 0;
    else
        temp_cone_weights = datarun.cone_weights(temp_cone_indices, RGC_indices(RGC));
    end
    
    % select a new location for the RGC
    temp_new_RGC_center = kosher_locations(RGC,:);
    % fine the nearby cones and assign weights;
    temp_distances = ipdm(temp_new_RGC_center, datarun.cone_centers);
    [sorted_cone_distances, sorted_indices] = sort(temp_distances, 'ascend');
    new_cone_indices = sorted_indices(1:temp_num_center_cones);
    
    % assign to new connectivity matrix and normalize
    new_connectivity_matrix(new_cone_indices, RGC) = temp_cone_weights ./ sum(temp_cone_weights);
    new_RGC_locations(RGC,:) = temp_new_RGC_center;
end

% remove the zero filled RGC weight vectors from the connectivity matrix.
% These zero filled weight vectors correspond to RGCs with a cone
% convergence less than min_cone_convergence or greater than
% max_cone_convergence
temp_keeper_indices = find(max(new_connectivity_matrix, [], 1) > 0);
new_connectivity_matrix = new_connectivity_matrix(:, temp_keeper_indices);

% visualization for check
RGC = 2;
figure(1)
clf
hold on
set(gcf, 'color', [0.2 0.2 0.2])
axis off    
%plot RGC centers
plot(new_RGC_locations(RGC,1), new_RGC_locations(RGC,2), 'wo')
% plot connection lines
temp_cone_indices = find(new_connectivity_matrix(:,RGC) > 0);
temp_num_cones = length(temp_cone_indices);
line_weights = new_connectivity_matrix(temp_cone_indices,RGC);
for cone = 1:temp_num_cones
    plot([new_RGC_locations(RGC,1) datarun.cone_centers(temp_cone_indices(cone),1)], [new_RGC_locations(RGC,2) datarun.cone_centers(temp_cone_indices(cone),2)], 'w', 'LineWidth', (abs(20*line_weights(cone))))
end
% plot cones
for cone = 1:num_cones
    if datarun.cone_types(cone) == 'L'
        plot(datarun.cone_centers(cone,1), datarun.cone_centers(cone,2), 'r.','Markersize', 10)
    elseif datarun.cone_types(cone) == 'M'
        plot(datarun.cone_centers(cone,1), datarun.cone_centers(cone,2), 'g.','Markersize', 10)
    elseif datarun.cone_types(cone) == 'S'
        plot(datarun.cone_centers(cone,1), datarun.cone_centers(cone,2), 'b.','Markersize', 10)
    else
        plot(datarun.cone_centers(cone,1), datarun.cone_centers(cone,2), 'k.','Markersize', 10)
    end
end
%title('cone mosaic','color', 'w', 'FontSize', 14, 'FontName', 'Helvetica')
%print(1, '/snle/home/gfield/Desktop/simulated_rfs','-dpdf')   

% extract convergence
convergences = zeros(1, num_RGCs-excluded_cell_counter);
for RGC = 1:(num_RGCs-excluded_cell_counter)
    convergences(RGC) = length(find(new_connectivity_matrix(:,RGC) > 0));
end
mean(convergences)
std(convergences)

% calculate the OI
temp_indices = compute_opponency_index(new_connectivity_matrix, datarun.cone_types);
current_std = std(temp_indices)



