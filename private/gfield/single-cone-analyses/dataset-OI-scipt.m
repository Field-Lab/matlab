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
hist(convergences,[0:0.5:50])

% calculate the OI
temp_indices = compute_opponency_index(new_connectivity_matrix, datarun.cone_types);
current_std = std(temp_indices)

% calculate the OI with random permutations
num_iters = 100;
temp_current_stds = zeros(1, num_iters);
for iters = 1:num_iters
    indices = randperm(num_cones);
    temp_indices = compute_opponency_index(new_connectivity_matrix, datarun.cone_types(indices));
    temp_current_stds(iters) = std(temp_indices);
end

% extract stats on permuted cone mosaics
mean(temp_current_stds)
std(temp_current_stds)
hist(temp_current_stds, [0.2:0.01:0.6])



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
