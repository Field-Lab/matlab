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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RGC locations
RGC_indices = find(datarun.cell_types == cell_type);
RGC_locations = datarun.rgc_COMs(RGC_indices, :);
%RGC_weights = datarun.cone_weights(:,RGC_indices);
num_RGCs = length(RGC_indices);
simulated_mosaic = false;



num_cycles = 20;

cycle_OI_SDs = zeros(1,num_cycles);
for cycle = 1:num_cycles

    % Plantain
    [RGC_rf_map, RGC_rf_weights, extras] = cone_to_RGC_connectivity(cone_locations, RGC_locations,...
        'center_amplitude', 10,'center_space_constant', 8, 'center_weight_sd_factor', 0.8,...
        'sharing_probability', 0.2, 'seed', 1111+num_cycles);

    % % Mango
    % [RGC_rf_map, RGC_rf_weights, extras] = cone_to_RGC_connectivity(cone_locations, RGC_locations,...
    %     'center_amplitude', 20,'center_space_constant', 16, 'center_weight_sd_factor', 0.4,...
    %     'sharing_probability', 0.2, 'seed', 6666);
    % 
    % % Plum
    % [RGC_rf_map, RGC_rf_weights, extras] = cone_to_RGC_connectivity(cone_locations, RGC_locations,...
    %     'center_amplitude', 12,'center_space_constant', 8, 'center_weight_sd_factor', 0.4,...
    %     'sharing_probability', 0.3, 'seed', 6666);
    % 
    % 
    % % blueberry
    % [RGC_rf_map, RGC_rf_weights, extras] = cone_to_RGC_connectivity(cone_locations, RGC_locations,...
    %     'center_amplitude', 12,'center_space_constant', 10, 'center_weight_sd_factor', 1.0,...
    %     'sharing_probability', 0.25, 'seed', 6666);
    % 
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

    length(RGC_locations(:,1));
    keepers;
    mean_sim_convergences = mean(sim_convergences);
    std_sim_convergences = std(sim_convergences);

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
    current_std = std(temp_indices);
    cycle_OI_SDs(cycle) = current_std;
end

figure(2)
mean(cycle_OI_SDs)
std(cycle_OI_SDs)
hist(cycle_OI_SDs)

