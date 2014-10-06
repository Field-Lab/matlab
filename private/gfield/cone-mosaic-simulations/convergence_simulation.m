%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% script for generating cone mosaics sampled by RGC mosaics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear

% generate mosaic of cones
sim_info.cones.l_proportion = 0.59;
sim_info.cones.mosaic_seed = 1111;
sim_info.cones.mosaic_method = 'exclusion_zone';  % or 'hexagonal'

% hexagonal cone mosaic
if strcmp(sim_info.cones.mosaic_method, 'hexagonal')
    field_buffer = 10; % imposes a boarder around the cone locations
    sim_info.cones.number = 60^2; % the squareroot must be an integer
    sim_info.cones.space_factor = 1.0;
    sim_info.cones.mosaic_SD = 0.2;
    temp_cone_locations = hexagonal_mosaic(sim_info.cones.number, 'spacing', sim_info.cones.space_factor,...
        'jitter_method', 'bounded_uniform', 'BU_scale', 0.3,'seed', sim_info.cones.mosaic_seed );
    % ensure all cone lications are positive for compatibility with data functions
    minimum_coordinate = min(min(temp_cone_locations));
    maximum_coordinate = max(max(temp_cone_locations));
    temp_cone_locations = temp_cone_locations + repmat((minimum_coordinate + field_buffer),sim_info.cones.number, 2);
    simrun.cones.centers = temp_cone_locations;
    simrun.stimulus.field_height = round(maximum_coordinate + minimum_coordinate + (2*field_buffer));
    simrun.stimulus.field_width = round(maximum_coordinate + minimum_coordinate + (2*field_buffer));
end

% exclusion cone mosaic
if strcmp(sim_info.cones.mosaic_method, 'exclusion_zone')
    field_buffer = 10;
    sim_info.cones.number = 50^2; % the squareroot must be an integer
    sim_info.cones.grid_size = [52 52];
    sim_info.cones.exclusion_mean = sqrt(3)/2;
    sim_info.cones.exclusion_sigma = exclusion_mean * 0.1;
    temp_cone_locations = make_serial_exclusion_mosaic(sim_info.cones.grid_size, sim_info.cones.number,...
                                                        sim_info.cones.exclusion_mean, sim_info.cones.exclusion_sigma,...
                                                        'verbose', true, 'seed', sim_info.cones.mosaic_seed);
    temp_cone_locations = temp_cone_locations + repmat(field_buffer, sim_info.cones.number, 2);
    simrun.cones.centers = temp_cone_locations;
    simrun.stimulus.field_height = sim_info.cones.grid_size(1) + (2*field_buffer);
    simrun.stimulus.field_width = simrun.stimulus.field_height;
end

%random
temp_types = assign_cone_types(simrun.cones.centers, sim_info.cones.l_proportion);
simrun.cones.types = temp_types;

% clumped
temp_types = assign_cone_types(simrun.cones.centers, sim_info.cones.l_proportion, 'clumping', true,...
                                'filter_sigma', 3.0, 'delta_pl', 0.9);
simrun.cones.types = temp_types;


plot_cone_mosaic(simrun)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% artifactual clumption
clumped_fraction = 0.1;
temp_simrun = make_artifactual_clumping(simrun, clumped_fraction);
% check num proportions after clumping         
test_proportion = length(find(cone_types == 0)) ./ num_cones
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

num_bins = 1200;
bin_width = 0.02;
l_cone_ids = find(cone_types == 0);
m_cone_ids = find(cone_types == 1);
% l cone analysis
% l cone to all cones
[full_cone_drp, bin_centers] = density_recovery_profile(simrun.cones.centers, num_bins, bin_width);
% l cones to l cones
[l_to_l_cone_drp, bin_centers] = density_recovery_profile(simrun.cones.centers(l_cone_ids,:), num_bins, bin_width);
% l cones to m cones
[l_to_m_cone_drp, bin_centers] = density_recovery_profile(simrun.cones.centers(m_cone_ids,:), num_bins, bin_width,...
    'reference_centers', simrun.cones.centers(l_cone_ids,:));

figure(3)
plot(bin_centers, full_cone_drp, 'k')
hold on
plot(bin_centers, l_to_l_cone_drp, 'r')
plot(bin_centers, l_to_m_cone_drp, 'g')

figure(11)
clf
hold on
plot(bin_centers, l_to_l_cone_drp ./ sum(l_to_l_cone_drp), 'r')
plot(bin_centers, l_to_m_cone_drp ./ sum(l_to_m_cone_drp) , 'g')
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cone packing factor
% temp_areas = pi*bin_width^2*(2*(1:num_bins)-1);
% expected_counts = num_cones * temp_areas .* cone_density;
% observed_counts = full_cone_drp * num_cones .* temp_areas;
% cutoff_index = find(observed_counts > expected_counts, 1);
% 
% dip_volume = sum(expected_counts(1:cutoff_index-1) - observed_counts(1:cutoff_index-1)) ./ num_cones;
% 
% exclusion_radius_squared = dip_volume ./ (pi * cone_density);
% max_radius_squared = sqrt(4/3) / cone_density;
% 
% packing_factor = exclusion_radius_squared / max_radius_squared
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate a mosaic of RGCs
sim_info.RGCs.mean_convergence = 11;
sim_info.RGCs.mosaic_seed = 2222;
sim_info.RGCs.mosaic_method = 'exclusion_zone';

% hexagonal RGC mosaic
if strcmp(sim_info.RGCs.mosaic_method, 'hexagonal')
    sim_info.RGCs.number = (round(sqrt(sim_info.cones.number ./ sim_info.RGCs.mean_convergence)))^2;
    sim_info.RGCs.mosaic_SD = 0.65;
    sim_info.RGCs.space_factor = sim_info.cones.space_factor * sqrt(sim_info.cones.number ./ sim_info.RGCs.number);
    temp_RGC_locations = hexagonal_mosaic(sim_info.RGCs.number, 'spacing', sim_info.RGCs.space_factor,...
        'jitter_method', 'bounded_uniform', 'BU_scale', 2.0, 'seed', sim_info.RGCs.mosaic_seed);
    simrun.rgc_COMs = temp_RGC_locations;
end

% exclusion zone RGC mosaic
if strcmp(mosaic_method, 'exclusion_zone')
    num_RGCs = (round(sqrt(num_cones ./ mean_convergence)))^2
    RGC_grid_size = cone_grid_size;
    % to calc the exclusion radius start with that for a hex lattice and
    % scale down by some fraction, i.e. 0.9.
    exclusion_mean = sqrt(sqrt(4/3) ./ (num_RGCs ./ prod(RGC_grid_size))) * 0.9;
    exclusion_sigma = exclusion_mean * 0.1;
    verbose = true;
    RGC_locations = make_serial_exclusion_mosaic(RGC_grid_size, num_RGCs, exclusion_mean, exclusion_sigma, 'verbose', verbose, 'seed', RGC_mosaic_seed);
end

% plot cone mosiac and RGC mosaic
figure(1)
clf
hold on
for cone = 1:num_cones
    if cone_types(cone) == 0
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% loop through the the responses

% seed the random number generator
rand('twister', 2222);
cone_resample_probability = 0.005;
miss_cone_probability = 0.3;
% initialze RGC_to_cone_map
RGCs_to_cones_map = zeros(num_RGCs);
assigned_cones = zeros(num_cones, 1);
cone_counter = 1; 
diverging_cones = 0; % initialize a counter that keeps track of the number of cones that diverge to >1 RGC
while cone_counter <= num_cones
    RGC_index = ceil(num_RGCs * rand(1));  % randomly select an RGC
    %fprintf('RGC # %d \n', RGC_index)
    if mod(cone_counter, 1000) == 0
        fprintf('%d unique cones assign to RGCs \n', cone_counter)
    end
    temp_RGC_location = RGC_locations(RGC_index, :);    % get the location of this RGC
    temp_RGC_location = repmat(temp_RGC_location, num_cones, 1);    % bookkeeping
    temp_distances = sqrt((temp_RGC_location(:,1) - cone_locations(:,1)).^2 + (temp_RGC_location(:,2) - cone_locations(:,2)).^2); % compute distances to each cone
    [ascending_distance, cone_list] = sort(temp_distances, 'ascend');   % ascenting sort of these distances
    free_cone = true;    % if true, then cone of interest hasn't been assigned to an RGC
    temp_index = 1; % index to progress through "cone_list:
    while free_cone
        check_cone_assignment = find(assigned_cones == cone_list(temp_index)); % has this cone been assigned to any RGC?
        if isempty(check_cone_assignment)    % if not, assign to RGC of interest
            unassigned_indices = find(RGCs_to_cones_map(RGC_index,:) == 0); % find the slots that haven't been assigned to cones
            RGCs_to_cones_map(RGC_index, unassigned_indices(1)) = cone_list(temp_index); % assign this cone to the first empty slow of this RGC
            assigned_cones(cone_counter) = cone_list(temp_index); % add this cone to the list of assigned cones
            cone_counter = cone_counter + 1; % increment the cones that have been assigned
            free_cone = false; % escape from while loop -- cone is no longer free, it has been assigned to an RGC
            %fprintf('cone # %d was assigned to RGC %d \n', cone_list(temp_index), RGC_index)
        else
            % has the cone been assigned to the RGC of interest?
            already_listed = find(RGCs_to_cones_map(RGC_index,:) == cone_list(temp_index)); 
            if ~isempty(already_listed)  % if it has been assigned to RGC of interest, go to next closest cone
                temp_index = temp_index + 1;
                %fprintf('incrementing to next cone \n')
            else    % cone has been assigned to a different RGC
                random_number = rand(1); % generate a random number to determine course of action
                if random_number <= cone_resample_probability
                    unassigned_indices = find(RGCs_to_cones_map(RGC_index,:) == 0); % find the slots that haven't been assigned to cones   
                    RGCs_to_cones_map(RGC_index, unassigned_indices(1)) = cone_list(temp_index); % assign this cone to the first empty slow of this RGC
                    free_cone = false; % escape from while loop -- cone has been assigned to more than one RGC
                    diverging_cones = diverging_cones + 1; % increment diverging_cones counter
                    %fprintf('cone # %d was assigned to RGC %d \n', cone_list(temp_index), RGC_index)
                    %fprintf('cone # %d is connected to more than one RGC \n', cone_list(temp_index))
                elseif (random_number >= cone_resample_probability) && (random_number <= (cone_resample_probability + miss_cone_probability))
                    %fprintf('RGC # %d did not connect to a cone on this iteration \n', RGC_index)
                    break % this RGC misses a chance to collect signals from a cone
                else    % move on to next cone
                    temp_index = temp_index + 1;
                    %fprintf('moving to next closest cone \n')
                end
            end
        end
    end
end

RGC_rf_map = cell(num_RGCs, 1); % initialize cell array
for RGC = 1:num_RGCs
    temp_rf = RGCs_to_cones_map(RGC,:);
    temp_zero_indices = find(temp_rf == 0);
    temp_rf = temp_rf(1:(temp_zero_indices(1)-1));
    RGC_rf_map{RGC} = temp_rf;
end


figure(1)
clf
hold on
%plot RGC centers
plot(RGC_locations(1,:), RGC_locations(2,:), 'w.')
% plot connection lines
for RGC = 1:num_RGCs
    temp_rf = RGC_rf_map{RGC};
    temp_num_cones = length(temp_rf);
    for cone = 1:temp_num_cones
        plot([RGC_locations(RGC,1) cone_locations(temp_rf(cone),1)], [RGC_locations(RGC,2) cone_locations(temp_rf(cone), 2)], 'w')
    end
end
% plot cones
for cone = 1:num_cones
    if cone_types(cone) == 0
        plot(cone_locations(cone,1), cone_locations(cone,2), 'r.')
    else
        plot(cone_locations(cone,1), cone_locations(cone,2), 'g.')
    end
end
set(gcf, 'color', [0.5 0.5 0.5])
%title('cone mosaic','color', 'w', 'FontSize', 14, 'FontName', 'Helvetica')
axis off    
%print(1, '/snle/home/gfield/Desktop/simulated_rfs','-dpdf')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% extract some statistics
convergences = zeros(num_RGCs, 1);
for RGC = 1:num_RGCs
    temp_cone_num = length(RGC_rf_map{RGC});
    convergences(RGC) = temp_cone_num;
end
figure(2)
hist(convergences,[0:1:30])
%print(2, '/snle/home/gfield/Desktop/convergence','-dpdf')
% divergence
cone_divergence = zeros(num_cones, 1);
for cone = 1:num_cones
    for RGC = 1:num_RGCs
        temp_rf = RGC_rf_map{RGC};
        temp_index = find(temp_rf == cone);
        if ~isempty(temp_index)
            cone_divergence(cone) = cone_divergence(cone) + 1;
        end
    end
end
figure(3)
hist(cone_divergence, [1:10])
%print(3, '/snle/home/gfield/Desktop/divergence','-dpdf')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simulate from the cones

[RGC_rf_map, RGC_rf_weights, extras] = cone_to_RGC_connectivity(cone_locations, RGC_locations, 'seed', 5555);

% define region of interest to avoid boarder effects.
buffer_size = 5*(sqrt(3)/2);
RGC_grid_size(1) = cone_grid_size(1) - buffer_size; % introduce a buffer of ~10 cones
RGC_grid_size(2) = cone_grid_size(2) - buffer_size;
x_indices = find(RGC_locations(:,1) > buffer_size & RGC_locations(:,1) < (cone_grid_size(1) - buffer_size));
y_indices = find(RGC_locations(:,2) > buffer_size & RGC_locations(:,2) < (cone_grid_size(1) - buffer_size));
roi_RGC_indices = intersect(x_indices, y_indices);

roi_RGC_locations = RGC_locations([roi_RGC_indices],:);
roi_RGC_rf_map = RGC_rf_map(roi_RGC_indices);
roi_RGC_rf_weights = RGC_rf_weights(roi_RGC_indices);
num_roi_RGCs = length(roi_RGC_indices);

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
    if cone_types(cone) == 0
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
convergences = zeros(num_roi_RGCs, 1);
for RGC = 1:num_roi_RGCs
    temp_cone_num = length(roi_RGC_rf_map{RGC});
    convergences(RGC) = temp_cone_num;
end
figure(2)
hist(convergences,[0:5:30])
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
    
%fprintf('there were %d cone that were not sampled by an RGC \n', unsampled_cones);    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pack simulation into a datarun for export

clear save_dir datarun_sim datarun

connectivity_matrix = zeros(num_cones, num_roi_RGCs);
for RGC = 1:num_roi_RGCs
    temp_rf = zeros(num_cones,1);
    temp_rf(roi_RGC_rf_map{RGC}) = roi_RGC_rf_weights{RGC};
    connectivity_matrix(:,RGC) = temp_rf;
end
%normalization
for RGC = 1:num_roi_RGCs
    center_cone_indices = find(connectivity_matrix(:,RGC) > 0);
    total_center_input = sum(connectivity_matrix(center_cone_indices,RGC));
    connectivity_matrix(:,RGC) = connectivity_matrix(:,RGC) ./ total_center_input;
end


clear new_cone_types
for cone = 1:num_cones
    if cone_types(cone) == 0
        new_cone_types(cone) = 'L';
    else
        new_cone_types(cone) = 'M';
    end
end

datarun_sim.cones.weights = connectivity_matrix;
datarun_sim.cones.centers = cone_locations(:,1:2);
datarun_sim.cones.types = new_cone_types;
datarun_sim.cones.rf_fits.center = roi_RGC_locations;


% test opponency index calculation
temp_indices = compute_opponency_index(connectivity_matrix, new_cone_types);
current_std = std(temp_indices)

bins = -1:0.1:1;
oi_hist = hist(temp_indices, bins);
bar(bins, oi_hist)

normalized_oi_width = mean(current_std ./ random_connections.sim_stds)
normalized_oi_width_error = std(current_std ./ random_connections.sim_stds)

% poor-mans generation of datarun structure from simulation
datarun.cone_centers = cone_locations;
datarun.rgc_COMs = roi_RGC_locations;
datarun.cone_types = new_cone_types;
datarun.cone_weights = connectivity_matrix;


