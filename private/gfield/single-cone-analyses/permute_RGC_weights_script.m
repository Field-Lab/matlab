%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% extract information from the data

datarun = import_single_cone_data([], 'plantain');
%datarun = import_single_cone_data([], 'mango');

cell_type = 4;
RGC_min_convergence = 5;
RGC_max_convergence = 30;
radius_scale = 2;

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



temp_RGC_indices = find(datarun.cell_types == cell_type);
num_RGCs = length(temp_RGC_indices);
num_cones = length(datarun.cone_centers);
temp_connect_matrix = datarun.cone_weights(:, temp_RGC_indices);

% get distances between cones and RGCs of interest
distance_matrix = ipdm(datarun.rgc_COMs(temp_RGC_indices,:), datarun.cone_centers);

connectivity_matrix = zeros(num_cones, num_RGCs);
for RGC = 1:num_RGCs
    temp_cone_indices = find(distance_matrix(RGC, :) < (radius_scale * datarun.rgc_fit_info(temp_RGC_indices(RGC), 3)));
    new_cone_indices = find(datarun.cone_weights(temp_cone_indices, temp_RGC_indices(RGC)) > 0);
    temp_cone_weights = temp_connect_matrix(temp_cone_indices(new_cone_indices), RGC);
    temp_norm_factor = sum(temp_cone_weights);
    temp_norm_weights = temp_cone_weights ./ temp_norm_factor;
    connectivity_matrix(temp_cone_indices(new_cone_indices), RGC) = temp_norm_weights;
end

% creat an ROI based on convergence
roi_RGC_indices = [];
for RGC = 1:num_RGCs
    if (length(find(connectivity_matrix(:,RGC))) >= RGC_min_convergence) && (length(find(connectivity_matrix(:,RGC))) <= RGC_max_convergence);
        roi_RGC_indices = [roi_RGC_indices, RGC];
    end
end
roi_RGC_COMs = datarun.rgc_COMs(temp_RGC_indices(roi_RGC_indices),:);
roi_connectivity_matrix = connectivity_matrix(:, roi_RGC_indices);
roi_RGC_locations = datarun.rgc_COMs(temp_RGC_indices(roi_RGC_indices),:);
roi_distance_matrix = distance_matrix(roi_RGC_indices,:);
num_roi_RGCs = length(roi_RGC_indices);


% calculate the OI
temp_indices = compute_opponency_index(roi_connectivity_matrix, datarun.cone_types);
data_OI_SD = std(temp_indices)

set(0,'DefaultAxesFontSize', 18, 'DefaultAxesFontName', 'Helvetica')

% plot distribution of OIs
[purity_hist, hist_bins] = hist(temp_indices, [-1.0:0.2:1.0]);
bar(hist_bins, purity_hist)
axis([-1.0 1.0 0 25])
print(1, '/snle/home/gfield/Desktop/off-par-hist','-dpdf')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% permute RGC Weights from one RGC to another

num_shuffles = 100;
for shuffle = 1:num_shuffles;
    %RGC_indices = 1:num_roi_RGCs;
    shuffled_RGC_indices = randperm(num_roi_RGCs);
    shuffled_connectivity_matrix = zeros(num_cones, num_roi_RGCs);
    for RGC = 1:num_roi_RGCs
        % get the weights for the RGC
        temp_cone_indices = find(roi_distance_matrix(RGC, :) < (radius_scale * datarun.rgc_fit_info(temp_RGC_indices(roi_RGC_indices(RGC)), 3)));
        temp_weights = roi_connectivity_matrix(temp_cone_indices, RGC);
        % order the weights according to distance from rgc COM
        temp_dists = roi_distance_matrix(RGC, temp_cone_indices);
        [sorted_dists, old_cone_indices] = sort(temp_dists, 'ascend');
        ordered_weights = temp_weights(old_cone_indices);
        % get a new RGC COM
        new_dists = roi_distance_matrix(shuffled_RGC_indices(RGC),:);
        [sorted_dists, new_cone_indices] = sort(new_dists, 'ascend');
        new_cone_indices = new_cone_indices(1:length(ordered_weights));
        % assign ordered weights to the closest cones about the new RGC location
        shuffled_connectivity_matrix(new_cone_indices, shuffled_RGC_indices(RGC)) = ordered_weights;
    end

    temp_indices = compute_opponency_index(shuffled_connectivity_matrix, datarun.cone_types);
    current_std = std(temp_indices);
    shuffled_OI_SDs(shuffle) = current_std;
end

mean(shuffled_OI_SDs)
std(shuffled_OI_SDs)
[permuted_RGC_hist, hist_bins] = hist(shuffled_OI_SDs, [0.1:0.01:0.6]);


% calculate the OI with random permutations
num_iters = 100;
temp_current_stds = zeros(1, num_iters);
for iters = 1:num_iters
    indices = randperm(num_cones);
    temp_indices = compute_opponency_index(roi_connectivity_matrix, datarun.cone_types(indices));
    temp_current_stds(iters) = std(temp_indices);
end

[permuted_cone_hist, hist_bins] = hist(temp_current_stds, hist_bins);

figure(1)
clf
hold on
bar(hist_bins, (-1*permuted_RGC_hist), 'm')
bar(hist_bins, permuted_cone_hist, 'b', 'k')
plot([data_OI_SD data_OI_SD], [-25 25], 'b')
axis([0.2 0.6001 -25 25])
print(1, '/snle/home/gfield/Desktop/growth-hist-2','-dpdf')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% shuffle RGCs to new (unconstrained) locations
% Each RGC samples from the cone mosaic independently of all other RGCs
jitter_scale = 4;
radius_scaler = 2;

temp_RGC_indices = find(datarun.cell_types == cell_type);

num_shuffles = 100;
free_RGC_shuffle_OI_SD = zeros(1,num_shuffles);
for shuffle = 1:num_shuffles

    % precompute the new RGC locations
    %rand('twister', 1111+shuffle);
    location_jitter = jitter_scale * rand(num_RGCs,2);
    roi_cones = find(datarun.cone_roi == 1);
    num_roi_cones = length(roi_cones);
    rand_cone_indices = randperm(num_roi_cones);
    kosher_locations = datarun.cone_centers(roi_cones(rand_cone_indices(1:num_RGCs)),:) + location_jitter;

    shuffled_connectivity_matrix = zeros(num_cones, num_RGCs);
    new_RGC_locations = zeros(num_RGCs,2);
    excluded_cell_counter = 0;
    for RGC = 1:length(temp_RGC_indices)

        % get number of center cones and associated weights
        temp_RGC_center = datarun.rgc_COMs(temp_RGC_indices(RGC),:);
        temp_RGC_CA = datarun.rgc_fit_info(temp_RGC_indices(RGC),3) * radius_scaler;
        temp_distances = ipdm(temp_RGC_center, datarun.cone_centers);
        close_cone_indices = find(temp_distances <= temp_RGC_CA);
        positive_cone_indices = find(datarun.cone_weights(:, temp_RGC_indices(RGC)) > 0);
        temp_cone_indices = intersect(close_cone_indices, positive_cone_indices);
        temp_num_center_cones = length(temp_cone_indices);

        clear temp_cone_weights
        if temp_num_center_cones <= RGC_min_convergence || temp_num_center_cones >= RGC_max_convergence
            excluded_cell_counter = excluded_cell_counter + 1;
            temp_cone_weights(1:temp_num_center_cones) = 0;
        else
            temp_cone_weights = datarun.cone_weights(temp_cone_indices, temp_RGC_indices(RGC));
        end

        % select a new location for the RGC
        temp_new_RGC_center = kosher_locations(RGC,:);
        % find the nearby cones and assign weights;
        temp_distances = ipdm(temp_new_RGC_center, datarun.cone_centers);
        [sorted_cone_distances, sorted_indices] = sort(temp_distances, 'ascend');
        new_cone_indices = sorted_indices(1:temp_num_center_cones);

        % assign to new connectivity matrix and normalize
        shuffled_connectivity_matrix(new_cone_indices, RGC) = temp_cone_weights ./ sum(temp_cone_weights);
        new_RGC_locations(RGC,:) = temp_new_RGC_center;
    end

    % remove the zero filled RGC weight vectors from the connectivity matrix.
    % These zero filled weight vectors correspond to RGCs with a cone
    % convergence less than min_cone_convergence or greater than
    % max_cone_convergence
    temp_keeper_indices = find(max(shuffled_connectivity_matrix, [], 1) > 0);
    shuffled_connectivity_matrix = shuffled_connectivity_matrix(:, temp_keeper_indices);

    % extract convergence
    convergences = zeros(1, num_RGCs - excluded_cell_counter);
    for RGC = 1:(num_RGCs - excluded_cell_counter)
        convergences(RGC) = length(find(shuffled_connectivity_matrix(:,RGC) > 0));
    end
    mean(convergences);
    std(convergences);

    % calculate the OI
    temp_indices = compute_opponency_index(shuffled_connectivity_matrix, datarun.cone_types);
    current_std = std(temp_indices);
    free_RGC_shuffle_OI_SD(shuffle) = current_std;
end

mean(free_RGC_shuffle_OI_SD)
std(free_RGC_shuffle_OI_SD)
hold on

[free_RGC_shuffle_hist, hist_bins] = hist(free_RGC_shuffle_OI_SD, hist_bins);
bar(hist_bins, free_RGC_shuffle_hist, 'c')

% visualization for check
% RGC = 2;
% figure(4)
% clf
% hold on
% set(gcf, 'color', [0.2 0.2 0.2])
% axis off    
% %plot RGC centers
% plot(new_RGC_locations(RGC,1), new_RGC_locations(RGC,2), 'wo')
% % plot connection lines
% temp_cone_indices = find(new_connectivity_matrix(:,RGC) > 0);
% temp_num_cones = length(temp_cone_indices);
% line_weights = new_connectivity_matrix(temp_cone_indices,RGC);
% for cone = 1:temp_num_cones
%     plot([new_RGC_locations(RGC,1) datarun.cone_centers(temp_cone_indices(cone),1)], [new_RGC_locations(RGC,2) datarun.cone_centers(temp_cone_indices(cone),2)], 'w', 'LineWidth', (abs(20*line_weights(cone))))
% end
% % plot cones
% for cone = 1:num_cones
%     if datarun.cone_types(cone) == 'L'
%         plot(datarun.cone_centers(cone,1), datarun.cone_centers(cone,2), 'r.','Markersize', 10)
%     elseif datarun.cone_types(cone) == 'M'
%         plot(datarun.cone_centers(cone,1), datarun.cone_centers(cone,2), 'g.','Markersize', 10)
%     elseif datarun.cone_types(cone) == 'S'
%         plot(datarun.cone_centers(cone,1), datarun.cone_centers(cone,2), 'b.','Markersize', 10)
%     else
%         plot(datarun.cone_centers(cone,1), datarun.cone_centers(cone,2), 'k.','Markersize', 10)
%     end
% end
%title('cone mosaic','color', 'w', 'FontSize', 14, 'FontName', 'Helvetica')
%print(1, '/snle/home/gfield/Desktop/simulated_rfs','-dpdf') 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reassining weights to same RGC
% shuffle RGCs to RGC locations
jitter_scale = 4;
radius_scaler = 2;

temp_RGC_indices = find(datarun.cell_types == cell_type);


shuffled_connectivity_matrix = zeros(num_cones, num_RGCs);
new_RGC_locations = zeros(num_RGCs,2);
excluded_cell_counter = 0;
for RGC = 1:length(temp_RGC_indices)

    % get number of center cones and associated weights
    temp_RGC_center = datarun.rgc_COMs(temp_RGC_indices(RGC),:);
    temp_RGC_CA = datarun.rgc_fit_info(temp_RGC_indices(RGC),3) * radius_scaler;
    temp_distances = ipdm(temp_RGC_center, datarun.cone_centers);
    close_cone_indices = find(temp_distances <= temp_RGC_CA);
    positive_cone_indices = find(datarun.cone_weights(:, temp_RGC_indices(RGC)) > 0);
    temp_cone_indices = intersect(close_cone_indices, positive_cone_indices);
    temp_num_center_cones = length(temp_cone_indices);

    clear temp_cone_weights
    if temp_num_center_cones <= RGC_min_convergence || temp_num_center_cones >= RGC_max_convergence
        excluded_cell_counter = excluded_cell_counter + 1;
        temp_cone_weights(1:temp_num_center_cones) = 0;
    else
        temp_cone_weights = datarun.cone_weights(temp_cone_indices, temp_RGC_indices(RGC));
    end


    % assign to new connectivity matrix and normalize
    new_cone_indices = close_cone_indices(1:temp_num_center_cones);
    shuffled_connectivity_matrix(new_cone_indices, RGC) = temp_cone_weights ./ sum(temp_cone_weights);
    new_RGC_locations(RGC,:) = temp_new_RGC_center;
end

% remove the zero filled RGC weight vectors from the connectivity matrix.
% These zero filled weight vectors correspond to RGCs with a cone
% convergence less than min_cone_convergence or greater than
% max_cone_convergence
temp_keeper_indices = find(max(shuffled_connectivity_matrix, [], 1) > 0);
shuffled_connectivity_matrix = shuffled_connectivity_matrix(:, temp_keeper_indices);

% extract convergence
convergences = zeros(1, num_RGCs - excluded_cell_counter);
for RGC = 1:(num_RGCs - excluded_cell_counter)
    convergences(RGC) = length(find(shuffled_connectivity_matrix(:,RGC) > 0));
end
mean(convergences);
std(convergences);

% calculate the OI
temp_indices = compute_opponency_index(shuffled_connectivity_matrix, datarun.cone_types);
reweighted_RGCs_OI_SD = std(temp_indices);


plot([reweighted_RGCs_OI_SD, reweighted_RGCs_OI_SD], [-30 30], 'm')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% permute cone weights with in the center.

num_shuffles = 100;
permuted_weights_OI_SDs = zeros(1,num_shuffles);
for shuffle = 1:num_shuffles;
    %RGC_indices = 1:num_roi_RGCs;
    shuffled_RGC_indices = randperm(num_roi_RGCs);
    shuffled_connectivity_matrix = zeros(num_cones, num_roi_RGCs);
    for RGC = 1:num_roi_RGCs
        % get the weights for the RGC
        temp_cone_indices = find(roi_distance_matrix(RGC, :) < (radius_scale * datarun.rgc_fit_info(temp_RGC_indices(roi_RGC_indices(RGC)), 3)));
        temp_weights = roi_connectivity_matrix(temp_cone_indices, RGC);
        new_weight_order = randperm(length(temp_cone_indices));
        
        shuffled_connectivity_matrix(temp_cone_indices, RGC) = temp_weights(new_weight_order);
    end

    temp_indices = compute_opponency_index(shuffled_connectivity_matrix, datarun.cone_types);
    current_std = std(temp_indices);
    permuted_weights_OI_SDs(shuffle) = current_std;
end

mean(permuted_weights_OI_SDs)
std(permuted_weights_OI_SDs)
[permuted_local_weights_hist, hist_bins] = hist(permuted_weights_OI_SDs, [0.1:0.01:0.6]);

bar(hist_bins, permuted_local_weights_hist, 'b')
    
    
    