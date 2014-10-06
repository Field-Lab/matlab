datarun = load_data('/snle/lab/Experiments/Array/Analysis/2008-08-27-5/data003/data003/data003')
datarun = load_params(datarun,struct('verbose',1));  
datarun = load_sta(datarun,'load_sta',[]);


datarun = import_single_cone_data(datarun, 'plantain');
verbose = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Acquire the distribution of cone weights
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

radius_scaler = 2.0;
cell_type = 4;
RGC_min_convergence = 5;
RGC_max_convergence = 25;

%cell_type_indices = find(datarun.cell_types == cell_type);
%length(cell_type_indices)

[connectivity, new_datarun, connectivity_extras] = extract_connectivity(datarun, cell_type, 'remove_cones', 'SU', 'required_sign', 'positive',...
                                            'min_convergence', RGC_min_convergence, 'max_convergence', RGC_max_convergence,...
                                            'max_radius', radius_scaler, 'normalize', true);


% Get cell_array of ordered cone weights
distance_matrix = ipdm(connectivity_extras.cone_locations, connectivity_extras.roi_RGC_locations);

[num_cones, num_RGCs] = size(connectivity);                                        

% generated an array of ordered connection distributions
ordered_weights = cell(1, RGC_max_convergence);
for RGC = 1:num_RGCs
    temp_indices = find(connectivity(:,RGC) > 0);
    temp_weights = connectivity(temp_indices, RGC);
    temp_dists = distance_matrix(temp_indices, RGC);
    [temp_sorted_dists, temp_sorted_indices] = sort(temp_dists, 'ascend');
    temp_weights = temp_weights(temp_sorted_indices);
    
    for ord = 1:length(temp_weights)
        temp_ord_weights = ordered_weights{ord};
        temp_ord_weights = [temp_ord_weights, temp_weights(ord)];
        ordered_weights{ord} = temp_ord_weights;
    end
end


% observe ordered weight fall off
figure(1)
clf
hold on
for ord = 1:length(ordered_weights)
    temp_weights = ordered_weights{ord};
    errorbar(ord, mean(temp_weights), std(temp_weights)./sqrt(length(temp_weights)))
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


                                        
all_cone_locations = datarun.cone_centers;
all_cone_types = datarun.cone_types;

l_indices = find(all_cone_types == 'L');
m_indices = find(all_cone_types == 'M');

l_m_indices = sort([l_indices; m_indices]);

cone_locations = connectivity_extras.cone_locations;

if verbose
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
    %axis off
    %set(gcf, 'color', 'k')
    axis square
    hold off
end

cd /snle/lab/Experiments/Array/Shared/one/
mkdir rgc-growth
cd rgc-growth/

cone_types = connectivity_extras.cone_types;
clear RGC_growth_sd
remove_fields{1} = 'cell_ids';
remove_fields{2} = 'cell_types';
remove_fields{3} = 'cells_in_roi';
remove_fields{4} = 'cone_weights';
remove_fields{5} = 'rgc_COMs';
remove_fields{6} = 'rgc_fit_info';
num_iters = 100;
for iter = 1:num_iters
    tic
    % generate a mosaic of RGCs
    mean_convergence = 11;
    RGC_mosaic_seed = 2222+iter;
    num_RGCs = (round(sqrt((num_cones) ./ mean_convergence) * 1.0)).^2;
    RGC_grid_size = [225 150];
    % to calc the exclusion radius start with that for a hex lattice and
    % scale down by some fraction, i.e. 0.9.
    exclusion_mean = sqrt(sqrt(4/3) ./ (num_RGCs ./ prod(RGC_grid_size))) * 0.9;
    exclusion_sigma = exclusion_mean * 0.1;
    %RGC_locations = make_serial_exclusion_mosaic(RGC_grid_size, num_RGCs, exclusion_mean, exclusion_sigma, 'verbose', true, 'seed', RGC_mosaic_seed);

    
    %RGC_locations_matrix(iter,:,:) = RGC_locations;
    RGC_locations = squeeze(RGC_locations_matrix(iter,:,:));
    
    % shift RGC_locations
    RGC_locations(:,1) = RGC_locations(:,1) + 0;
    RGC_locations(:,2) = RGC_locations(:,2) + 100;

    % plot cone mosiac and RGC mosaic
    if verbose
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
        title('cone mosaic','color', 'w', 'FontSize', 14, 'FontName', 'Helvetica')
        axis off
        set(gcf, 'color', 'k')
        plot(RGC_locations(:,1), RGC_locations(:,2), 'wo')
        hold off
    end

    [RGC_rf_map, RGC_rf_weights, extras] = cone_to_RGC_connectivity(cone_locations, RGC_locations,...
    'center_amplitude', 7,'center_space_constant', 7, 'center_weight_sd_factor', 0.8,...
    'sharing_probability', 0.20, 'seed', 11110+iter);

    
    roi_RGC_indices = [];
    for RGC = 1:num_RGCs
        if (length(RGC_rf_map{RGC}) >= RGC_min_convergence) && (length(RGC_rf_map{RGC}) <= RGC_max_convergence);
            roi_RGC_indices = [roi_RGC_indices, RGC];
        end
    end


    roi_RGC_locations = RGC_locations([roi_RGC_indices],:);
    roi_RGC_rf_map = RGC_rf_map(roi_RGC_indices);
    roi_RGC_rf_weights = RGC_rf_weights(roi_RGC_indices);
    num_roi_RGCs = length(roi_RGC_indices)

    roi_connectivity_matrix = zeros(num_cones, num_roi_RGCs);
    for RGC = 1:num_roi_RGCs
        temp_rf = zeros(num_cones,1);
        temp_rf(roi_RGC_rf_map{RGC}) = roi_RGC_rf_weights{RGC};
        roi_connectivity_matrix(:,RGC) = temp_rf;
    end
    for RGC = 1:num_roi_RGCs
        center_cone_indices = find(roi_connectivity_matrix(:,RGC) > 0);
        total_center_input = sum(roi_connectivity_matrix(center_cone_indices,RGC));
        roi_connectivity_matrix(:,RGC) = roi_connectivity_matrix(:,RGC) ./ total_center_input;
    end

    % replace weights in simulated connectivity matrix with ordered weights
    % from data
    new_connectivity_matrix = zeros(num_cones, num_roi_RGCs);
    sim_distance_matrix = ipdm(cone_locations, roi_RGC_locations);
    for RGC = 1:num_roi_RGCs
        temp_indices = find(roi_connectivity_matrix(:,RGC));
        temp_distances = sim_distance_matrix(temp_indices, RGC);
        [sorted_temp_dists, sorted_indices] = sort(temp_distances, 'ascend');
        clear new_weights
        for cone = 1:length(temp_indices)
            temp_weights = ordered_weights{cone};
            if isempty(temp_weights)
                temp_weights = ordered_weights{18};
            end
            rand_index = floor(length(temp_weights)*rand(1))+1;
            new_weights(cone) = temp_weights(rand_index);
        end
        new_connectivity_matrix(temp_indices(sorted_indices),RGC) = new_weights;
    end


    % test opponency index calculation
    temp_indices = compute_opponency_index(new_connectivity_matrix, connectivity_extras.cone_types);
    current_std = std(temp_indices)

    bins = -1:0.1:1;
    oi_hist = hist(temp_indices, bins);
    bar(bins, oi_hist)

    RGC_growth_sd(iter) = current_std;
        
    temp_num_RGCs = length(RGC_locations_matrix(iter,:,1));
    temp_datarun = datarun;
    temp_datarun = rmfield(temp_datarun, remove_fields);
    temp_datarun.cell_ids(1:temp_num_RGCs) = 0;
    temp_datarun.cell_types(1:temp_num_RGCs) = cell_type;
    temp_datarun.cells_in_roi(1:temp_num_RGCs) = 0;
    temp_datarun.cone_weights = new_connectivity_matrix;
    temp_datarun.rgc_COMs = squeeze(RGC_locations_matrix(iter,:,:));
    temp_datarun.rgc_fit_info = zeros(temp_num_RGCs, 8);
    

    sim_datarun = temp_datarun;
    export_single_cone_data(sim_datarun,{1,2,3,4,5},[],['/snle/lab/Experiments/Array/Shared/one/rgc-growth' 'plantain' '/'])
    
    toc
end

[growth_hist, hist_bins] = hist(RGC_growth_sd, [0.2:0.01:0.6]);
bar(hist_bins, -1 * growth_hist, 'm')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check out distribution of connectivities
distance_list = [];
weight_list = [];
for RGC = 1:num_roi_RGCs
    temp_cone_indices = find(roi_connectivity_matrix(:,RGC) > 0);
    temp_distances = ipdm(roi_RGC_locations(RGC,:), cone_locations(temp_cone_indices,:));
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
plot(distance_list, weight_list, 'k.')
axis([0 20 0 0.7])

hist_bins = [0:1:19];
for bin = 1:length(hist_bins)
    temp_indices_one = find(distance_list < bin);
    temp_indices_two = find(distance_list >= bin-1);
    temp_indices = intersect(temp_indices_one, temp_indices_two);
    temp_weights = weight_list(temp_indices);
    binned_means(bin) = mean(temp_weights);
    binned_stds(bin) = std(temp_weights);
end

errorbar(hist_bins, binned_means, binned_stds)
%%%%%%%%%%%%%%%%%%%%%%%%%%%


figure(1)
clf
hold on
%plot RGC centers
plot(roi_RGC_locations(:,1), roi_RGC_locations(:,2), 'wo')
% plot connection lines
for RGC = 1:num_roi_RGCs
    temp_rf = roi_RGC_rf_map{RGC};
    temp_num_cones = length(temp_rf);
    line_weights = roi_RGC_rf_weights{RGC};
    for cone = 1:temp_num_cones
        plot([roi_RGC_locations(RGC,1) cone_locations(temp_rf(cone),1)], [roi_RGC_locations(RGC,2) cone_locations(temp_rf(cone), 2)], 'w', 'LineWidth', (abs(0.5*line_weights(cone))))
    end
end
% plot cones
for cone = 1:num_cones
    if cone_types(cone) == 'L'
        plot(cone_locations(cone,1), cone_locations(cone,2), 'r.','Markersize', 10)
    else
        plot(cone_locations(cone,1), cone_locations(cone,2), 'g.','Markersize', 10)
    end
end
set(gcf, 'color', [0.2 0.2 0.2])
%title('cone mosaic','color', 'w', 'FontSize', 14, 'FontName', 'Helvetica')
axis off    
%print(1, '/snle/home/gfield/Desktop/simulated_rfs','-dpdf')   


%%% extract some statistics
% actual RGC convergence
convergences = zeros(length(connectivity(1,:)), 1);
for RGC = 1:length(connectivity(1,:))
    temp_cone_num = length(find(connectivity(:,RGC)));
    convergences(RGC) = temp_cone_num;
end
figure(1)
hist(convergences,[0:2:30])
title('convergence')
mean(convergences)
std(convergences)


%%% extract some statistics
convergences = zeros(num_roi_RGCs, 1);
for RGC = 1:num_roi_RGCs
    temp_cone_num = length(roi_RGC_rf_map{RGC});
    convergences(RGC) = temp_cone_num;
end
figure(2)
hist(convergences,[0:2:30])
title('convergence')
%print(2, '/snle/home/gfield/Desktop/convergence','-dpdf')
% divergence
cone_divergence = zeros(num_cones, 1);
for cone = 1:num_cones
    for RGC = 1:num_roi_RGCs
        temp_rf = roi_RGC_rf_map{RGC};
        temp_index = find(temp_rf == cone);
        if ~isempty(temp_index)
            cone_divergence(cone) = cone_divergence(cone) + 1;
        end
    end
end
figure(3)
hist(cone_divergence, [0:5])
title('divergence')
%print(3, '/~Desktop/divergence','-dpdf')
mean(convergences)
std(convergences)


% test opponency index calculation
temp_indices = compute_opponency_index(roi_connectivity_matrix, cone_types);
current_std = std(temp_indices)

bins = -1:0.1:1;
oi_hist = hist(temp_indices, bins);
bar(bins, oi_hist)

% extract distributions of cone weights






