% Data set of potential interest
% Plantain: ON and OFF midgets
% Mango: ON midgets
% Plum: OFF midgets
% Peach: ON-midgets
% Blueberry: ON midgets and maybe OFF midgets (both are very complete)

datarun = import_single_cone_data([], 'plantain');
datarun = import_single_cone_data([], 'blueberry');

cell_type =4;
RGC_min_convergence = 5;
RGC_max_convergence = 24;

temp_RGC_indices = find(datarun.cell_types == cell_type);
num_RGCs = length(temp_RGC_indices);
num_cones = length(datarun.cone_centers);

temp_connect_matrix = datarun.cone_weights(:, temp_RGC_indices);

% get distances between cones and RGCs of interest
distance_matrix = ipdm(datarun.rgc_COMs(temp_RGC_indices,:), datarun.cone_centers);

radius_scale = 2;

connectivity_matrix = zeros(num_cones, num_RGCs);
for RGC = 1:num_RGCs
    temp_cone_indices = find(distance_matrix(RGC, :) < (radius_scale * datarun.rgc_fit_info(temp_RGC_indices(RGC), 3)));
      %datarun.cone_weights(temp_cone_indices, temp_RGC_indices(RGC))
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
num_roi_RGCs = length(roi_RGC_indices)


% Check out distribution of connectivities
distance_list = [];
weight_list = [];
for RGC = 1:num_roi_RGCs
    temp_cone_indices = find(roi_connectivity_matrix(:,RGC) > 0);
    temp_distances = ipdm(datarun.rgc_COMs(temp_RGC_indices(roi_RGC_indices(RGC)),:), datarun.cone_centers(temp_cone_indices,:));
    temp_cone_weights = roi_connectivity_matrix(temp_cone_indices, RGC);

    %plot(temp_distances, temp_cone_weights, 'k.')
    
    %figure(1)
%     clf
%     hold on
%     plot(datarun.rgc_COMs(temp_RGC_indices(roi_RGC_indices(RGC)), 1), datarun.rgc_COMs(temp_RGC_indices(roi_RGC_indices(RGC)), 2), 'r.')
%     plot(datarun.cone_centers(temp_cone_indices,1), datarun.cone_centers(temp_cone_indices,2), 'k.')
%     hold off
    
    distance_list = [distance_list, temp_distances];
    weight_list = [weight_list, temp_cone_weights'];
end
figure(10)
plot(distance_list, weight_list, 'k.')

hist_bins = [0:0.5:19];
for bin = 1:(length(hist_bins)-1)
    temp_indices_one = find(distance_list < hist_bins(bin+1));
    temp_indices_two = find(distance_list >= hist_bins(bin));
    temp_indices = intersect(temp_indices_one, temp_indices_two);
    temp_weights = weight_list(temp_indices);
    binned_means(bin) = mean(temp_weights);
    binned_stds(bin) = std(temp_weights);
    binned_weights{bin} = temp_weights(:);
end
figure(11)
errorbar(hist_bins, binned_means, binned_stds)

% progressive binning
bin_cylces = 32;
for cycle = 1:bin_cycles
    plot

% calculate the OI
temp_indices = compute_opponency_index(roi_connectivity_matrix, datarun.cone_types);
current_std = std(temp_indices)

% calculate the OI with random permutations
num_iters = 100;
temp_current_stds = zeros(1, num_iters);
for iters = 1:num_iters
    indices = randperm(num_cones);
    temp_indices = compute_opponency_index(roi_connectivity_matrix, datarun.cone_types(indices));
    temp_current_stds(iters) = std(temp_indices);
end

hist(temp_current_stds, [0.1:0.005:0.4])

    
bins = -1:0.1:1;
oi_hist = hist(temp_indices, bins);
bar(bins, oi_hist)


% extract convergence
convergences = zeros(1, num_RGCs);
for RGC = 1:num_roi_RGCs
    convergences(RGC) = length(find(roi_connectivity_matrix(:,RGC) > 0));
end
mean(convergences)
std(convergences)


hist(convergences,[0:2:30])





