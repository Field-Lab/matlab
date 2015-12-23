data_fruit = 'plantain';
cell_type = 4;
radius_scaler = 2;
RGC_min_convergence = 5;
RGC_max_convergence = 24;


datarun = import_single_cone_data([], data_fruit);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extract L and M cone information from mosaics

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


cone_locations = datarun.cone_centers;
cone_types = datarun.cone_types;
num_cones = length(cone_types);

l_indices = find(cone_types == 'L');
m_indices = find(cone_types == 'M');

% observe cone mosaic
figure(1)
clf
hold on
for cone = 1:num_cones
    if cone_types(cone) == 'L'
        plot(cone_locations(cone,1), cone_locations(cone,2), 'r.')
    else
        plot(cone_locations(cone,1), cone_locations(cone,2), 'g.')
    end
end
title('cone mosaic','color', 'w', 'FontSize', 14, 'FontName', 'Helvetica');
set(gcf, 'color', [0.5 0.5 0.5])
axis square
axis off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RGC locations
RGC_indices = find(datarun.cell_types == cell_type);
RGC_locations = datarun.rgc_COMs(RGC_indices, :);
%RGC_weights = datarun.cone_weights(:,RGC_indices);
num_RGCs = length(RGC_indices);
simulated_mosaic = false;


% restrict connectivity to consider only center cones
new_connectivity_matrix = zeros(num_cones, num_RGCs);
excluded_cell_counter = 0;
keeper_cells = 0;
clear convergences;
for RGC = 1:num_RGCs
    % get number of center cones and associated weights
    temp_RGC_center = datarun.rgc_COMs(RGC_indices(RGC),:);
    if simulated_mosaic == false
        temp_RGC_CA = datarun.rgc_fit_info(RGC_indices(RGC),3) * radius_scaler; % collecting distance for center coness
        temp_distances = ipdm(temp_RGC_center, datarun.cone_centers);
        close_cone_indices = find(temp_distances <= temp_RGC_CA);
        positive_cone_indices = find(datarun.cone_weights(:, RGC_indices(RGC)) > 0);
        temp_cone_indices = intersect(close_cone_indices, positive_cone_indices);
        temp_num_center_cones = length(temp_cone_indices);
    else
        temp_cone_indices = find(datarun.cone_weights(:,RGC) > 0);
        temp_num_center_cones = length(temp_cone_indices);
    end
  
    clear temp_cone_weights
    if temp_num_center_cones < RGC_min_convergence || temp_num_center_cones > RGC_max_convergence
        excluded_cell_counter = excluded_cell_counter + 1;
        temp_cone_weights(1:temp_num_center_cones) = 0;
    else
        temp_cone_weights = datarun.cone_weights(temp_cone_indices, RGC_indices(RGC));
        keeper_cells = 1+keeper_cells;
        convergences(keeper_cells) = temp_num_center_cones;
    end
    
    % assign to new connectivity matrix and normalize
    new_connectivity_matrix(temp_cone_indices, RGC) = temp_cone_weights ./ sum(temp_cone_weights);
end
% remove the zero filled RGC weight vectors from the connectivity matrix.
% These zero filled weight vectors correspond to RGCs with a cone
% convergence less than min_cone_convergence or greater than
% max_cone_convergence
temp_keeper_indices = find(max(new_connectivity_matrix, [], 1) > 0);
new_connectivity_matrix = new_connectivity_matrix(:, temp_keeper_indices);
roi_RGC_locations = RGC_locations([temp_keeper_indices],:);

[num_roi_cones, num_roi_RGCs] = size(new_connectivity_matrix);

mean_convergence = mean(convergences)
stdv_convergence = std(convergences)
hist(convergences,[0:2:100])

% cone and RGCs
figure(3)
clf
hold on
set(gcf, 'color', [0.5 0.5 0.5])
axis square
axis off
plot(roi_RGC_locations(:,1), roi_RGC_locations(:,2), 'wo')
% plot connection lines
for RGC = 1:num_roi_RGCs
    temp_indices = find(new_connectivity_matrix(:,RGC) > 0);
    temp_weights = new_connectivity_matrix(temp_indices,RGC);
    temp_num_cones = length(temp_indices);
    for cone = 1:temp_num_cones
        plot([roi_RGC_locations(RGC,1) cone_locations(temp_indices(cone),1)],...
            [roi_RGC_locations(RGC,2) cone_locations(temp_indices(cone), 2)],...
            'w', 'LineWidth', (abs(30.0*temp_weights(cone))))
        %abs(10.0*temp_weights(cone))
        %pause
    end
end
for cone = 1:num_cones
    if cone_types(cone) == 'L'
        plot(cone_locations(cone,1), cone_locations(cone,2), 'r.')
    else
        plot(cone_locations(cone,1), cone_locations(cone,2), 'g.')
    end
end
title('cone mosaic','color', 'w', 'FontSize', 14, 'FontName', 'Helvetica');
hold off


%%%%% ATH stuff %%%%%%

tmp=datarun.rgc_COMs;
myCells=tmp(:,1)>120&tmp(:,1)<220&tmp(:,2)>160&tmp(:,2)<210;
tmp=datarun.cone_centers;
myCones=tmp(:,1)>120&tmp(:,1)<220&tmp(:,2)>160&tmp(:,2)<210;

% RGC locations
RGC_indices = find(myCells);
RGC_locations = datarun.rgc_COMs(RGC_indices, :);
%RGC_weights = datarun.cone_weights(:,RGC_indices);
num_RGCs = length(RGC_indices);
simulated_mosaic = false;

myWeights=datarun.cone_weights(myCones,myCells);
myLocations=datarun.cone_centers(myCones,:);
myCellTypes=datarun.cell_types(myCells);
myConeTypes=datarun.cone_types(myCones);

fid=fopen('/Users/alexth/Desktop/cones_plantain_subset_coneTypes.txt','w')
for i=1:size(myConeTypes,1)
        fprintf(fid,'%c\n',myConeTypes(i));
end
fclose(fid)

fid=fopen('/Users/alexth/Desktop/cones_plantain_subset_cellTypes.txt','w')
for i=1:size(myCellTypes,1)
        fprintf(fid,'%d\n',myCellTypes(i));
end
fclose(fid)



fid=fopen('/Users/alexth/Desktop/cones_plantain_subset.txt','w')
for i=1:size(myWeights,1)
    for j=1:size(myWeights,2)
        fprintf(fid,'%f\t',myWeights(i,j));
    end
    fprintf(fid,'\n');
end
fclose(fid)
fid=fopen('/Users/alexth/Desktop/cones_plantain_subset_locations.txt','w')
for i=1:size(myLocations,1)
    for j=1:size(myLocations,2)
        fprintf(fid,'%f\t',myLocations(i,j));
    end
    fprintf(fid,'\n');
end
fclose(fid)


% restrict connectivity to consider only center cones
new_connectivity_matrix = zeros(num_cones, num_RGCs);
excluded_cell_counter = 0;
keeper_cells = 0;
clear convergences;
for RGC = 1:num_RGCs
    % get number of center cones and associated weights
    temp_RGC_center = datarun.rgc_COMs(RGC_indices(RGC),:);
    if simulated_mosaic == false
        temp_RGC_CA = datarun.rgc_fit_info(RGC_indices(RGC),3) * radius_scaler; % collecting distance for center coness
        temp_distances = ipdm(temp_RGC_center, datarun.cone_centers);
        close_cone_indices = find(temp_distances <= temp_RGC_CA);
        positive_cone_indices = find(datarun.cone_weights(:, RGC_indices(RGC)) > 0);
        temp_cone_indices = intersect(close_cone_indices, positive_cone_indices);
        temp_num_center_cones = length(temp_cone_indices);
    else
        temp_cone_indices = find(datarun.cone_weights(:,RGC) > 0);
        temp_num_center_cones = length(temp_cone_indices);
    end
  
    clear temp_cone_weights
    if temp_num_center_cones < RGC_min_convergence || temp_num_center_cones > RGC_max_convergence
        excluded_cell_counter = excluded_cell_counter + 1;
        temp_cone_weights(1:temp_num_center_cones) = 0;
    else
        temp_cone_weights = datarun.cone_weights(temp_cone_indices, RGC_indices(RGC));
        keeper_cells = 1+keeper_cells;
        convergences(keeper_cells) = temp_num_center_cones;
    end
    
    % assign to new connectivity matrix and normalize
    new_connectivity_matrix(temp_cone_indices, RGC) = temp_cone_weights ./ sum(temp_cone_weights);
end
% remove the zero filled RGC weight vectors from the connectivity matrix.
% These zero filled weight vectors correspond to RGCs with a cone
% convergence less than min_cone_convergence or greater than
% max_cone_convergence
temp_keeper_indices = find(max(new_connectivity_matrix, [], 1) > 0);
new_connectivity_matrix = new_connectivity_matrix(:, temp_keeper_indices);
roi_RGC_locations = RGC_locations([temp_keeper_indices],:);

[num_roi_cones, num_roi_RGCs] = size(new_connectivity_matrix);

mean_convergence = mean(convergences)
stdv_convergence = std(convergences)
hist(convergences,[0:2:100])

% cone and RGCs
figure(4)
clf
hold on
set(gcf, 'color', [0.5 0.5 0.5])
axis square
axis off
plot(roi_RGC_locations(:,1), roi_RGC_locations(:,2), 'wo')
% plot connection lines
for RGC = 1:num_roi_RGCs
    temp_indices = find(new_connectivity_matrix(:,RGC) > 0);
    temp_weights = new_connectivity_matrix(temp_indices,RGC);
    temp_num_cones = length(temp_indices);
    for cone = 1:temp_num_cones
        plot([roi_RGC_locations(RGC,1) cone_locations(temp_indices(cone),1)],...
            [roi_RGC_locations(RGC,2) cone_locations(temp_indices(cone), 2)],...
            'w', 'LineWidth', (abs(30.0*temp_weights(cone))))
        %abs(10.0*temp_weights(cone))
        %pause
    end
end
for cone = 1:num_cones
    if cone_types(cone) == 'L'
        plot(cone_locations(cone,1), cone_locations(cone,2), 'r.')
    else
        plot(cone_locations(cone,1), cone_locations(cone,2), 'g.')
    end
end
title('cone mosaic','color', 'w', 'FontSize', 14, 'FontName', 'Helvetica');
hold off





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now grow connections between RGC locations and cones

% Plantain
[RGC_rf_map, RGC_rf_weights, extras] = cone_to_RGC_connectivity(cone_locations, RGC_locations,...
    'center_amplitude', 10,'center_space_constant', 8, 'center_weight_sd_factor', 0.8,...
    'sharing_probability', 0.2, 'seed', 6666);

% Mango
[RGC_rf_map, RGC_rf_weights, extras] = cone_to_RGC_connectivity(cone_locations, RGC_locations,...
    'center_amplitude', 20,'center_space_constant', 16, 'center_weight_sd_factor', 0.4,...
    'sharing_probability', 0.2, 'seed', 6666);

% Plum
[RGC_rf_map, RGC_rf_weights, extras] = cone_to_RGC_connectivity(cone_locations, RGC_locations,...
    'center_amplitude', 12,'center_space_constant', 8, 'center_weight_sd_factor', 0.4,...
    'sharing_probability', 0.3, 'seed', 6666);


% blueberry
[RGC_rf_map, RGC_rf_weights, extras] = cone_to_RGC_connectivity(cone_locations, RGC_locations,...
    'center_amplitude', 12,'center_space_constant', 10, 'center_weight_sd_factor', 1.0,...
    'sharing_probability', 0.25, 'seed', 6666);

roi_RGC_indices = [];
keepers = 0;
clear sim_convergences
for RGC = 1:length(RGC_locations(:,1))
    if (length(RGC_rf_map{RGC}) >= RGC_min_convergence) && (length(RGC_rf_map{RGC}) <= RGC_max_convergence);
        roi_RGC_indices = [roi_RGC_indices, RGC];
        keepers = keepers +1;
        sim_convergences(keepers) = length(RGC_rf_map{RGC});
    end
end

length(RGC_locations(:,1))
keepers
mean_sim_convergences = mean(sim_convergences)
std_sim_convergences = std(sim_convergences)

roi_RGC_locations = RGC_locations([roi_RGC_indices],:);
roi_RGC_rf_map = RGC_rf_map(roi_RGC_indices);
roi_RGC_rf_weights = RGC_rf_weights(roi_RGC_indices);
num_roi_RGCs = length(roi_RGC_indices);

% generate a connectivity matrix
new_connectivity_matrix = zeros(num_cones, num_roi_RGCs);
for RGC = 1:num_roi_RGCs
    temp_rf = zeros(num_cones,1);
    temp_rf(roi_RGC_rf_map{RGC}) = roi_RGC_rf_weights{RGC};
    new_connectivity_matrix(:,RGC) = temp_rf;
end
%normalize weights
for RGC = 1:num_roi_RGCs
    center_cone_indices = find(new_connectivity_matrix(:,RGC) > 0);
    total_center_input = sum(new_connectivity_matrix(center_cone_indices,RGC));
    new_connectivity_matrix(:,RGC) = new_connectivity_matrix(:,RGC) ./ total_center_input;
end

% compute opponency indices
temp_indices = compute_opponency_index(new_connectivity_matrix, cone_types);
current_std = std(temp_indices)


% Check out distribution of connectivities
distance_list = [];
weight_list = [];
for RGC = 1:num_roi_RGCs
    temp_cone_indices = find(new_connectivity_matrix(:,RGC) > 0);
    temp_distances = ipdm(RGC_locations(roi_RGC_indices(RGC),:), datarun.cone_centers(temp_cone_indices,:));
    temp_cone_weights = new_connectivity_matrix(temp_cone_indices, RGC);

%     plot(temp_distances, temp_cone_weights, 'k.')
%     
%     figure(1)
%     clf
%     hold on
%     plot(RGC_locations(roi_RGC_indices(RGC), 1), RGC_locations(roi_RGC_indices(RGC), 2), 'r.')
%     plot(datarun.cone_centers(temp_cone_indices,1), datarun.cone_centers(temp_cone_indices,2), 'k.')
%     hold off
% 
    
    distance_list = [distance_list, temp_distances];
    weight_list = [weight_list, temp_cone_weights'];
end
figure(15)
plot(distance_list, weight_list, 'k.')

hist_bins = [0:0.5:19];
for bin = 1:(length(hist_bins)-1)
    temp_indices_one = find(distance_list < hist_bins(bin+1));
    temp_indices_two = find(distance_list >= hist_bins(bin));
    temp_indices = intersect(temp_indices_one, temp_indices_two);
    temp_weights = weight_list(temp_indices);
    binned_means(bin) = mean(temp_weights);
    binned_stds(bin) = std(temp_weights);
end
figure(16)
errorbar(hist_bins, binned_means, binned_stds)


% cone and RGCs
figure(8)
clf
hold on
set(gcf, 'color', [0.5 0.5 0.5])
axis square
axis off
plot(roi_RGC_locations(:,1), roi_RGC_locations(:,2), 'wo')
% plot connection lines
for RGC = 1:num_roi_RGCs
    temp_indices = find(new_connectivity_matrix(:,RGC) > 0);
    temp_weights = new_connectivity_matrix(temp_indices,RGC);
    temp_num_cones = length(temp_indices);
    for cone = 1:temp_num_cones
        plot([roi_RGC_locations(RGC,1) cone_locations(temp_indices(cone),1)],...
            [roi_RGC_locations(RGC,2) cone_locations(temp_indices(cone), 2)],...
            'w', 'LineWidth', (abs(30.0*temp_weights(cone))))
    end
end
for cone = 1:num_cones
    if cone_types(cone) == 'L'
        plot(cone_locations(cone,1), cone_locations(cone,2), 'r.')
    else
        plot(cone_locations(cone,1), cone_locations(cone,2), 'g.')
    end
end

title('cone mosaic','color', 'w', 'FontSize', 14, 'FontName', 'Helvetica');
hold off








