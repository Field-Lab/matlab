%%% To Dos:  Check normalization %%%


%%%% define cone locations %%%%
% generate mosaic of cones
l_proportion = 0.59;
cone_mosaic_seed = 1111;
mosaic_method = 'exclusion_zone';  % or 'hexagonal'

% exclusion cone mosaic
num_cones = 50^2; % the squareroot must be an integer
if strcmp(mosaic_method, 'exclusion_zone')
    cone_grid_size = [52 52];
    exclusion_mean = sqrt(3)/2;
    exclusion_sigma = exclusion_mean * 0.1;
    verbose = true;
    cone_locations = make_serial_exclusion_mosaic(cone_grid_size, num_cones, exclusion_mean, exclusion_sigma, 'verbose', verbose, 'seed', cone_mosaic_seed);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pL = 0.59;
filter_sigma = 3;
num_cycles = 100;
deltapL = [0.0];
mean_convergence = [11];
normalization = true;

mean_oi = zeros(length(mean_convergence), length(deltapL));
error_oi = zeros(length(mean_convergence), length(deltapL));

for converg_index = 1:length(mean_convergence)

    % exclusion zone RGC mosaic
    RGC_mosaic_seed = 2222;
    num_RGCs = (round(sqrt(num_cones ./ mean_convergence(converg_index))))^2;
    RGC_grid_size = cone_grid_size;
    % to calc the exclusion radius start with that for a hex lattice and
    % scale down by some fraction, i.e. 0.9.
    exclusion_mean = sqrt(sqrt(4/3) ./ (num_RGCs ./ prod(RGC_grid_size))) * 0.9;
    exclusion_sigma = exclusion_mean * 0.1;
    verbose = true;
    RGC_locations = make_serial_exclusion_mosaic(RGC_grid_size, num_RGCs, exclusion_mean, exclusion_sigma, 'verbose', verbose, 'seed', RGC_mosaic_seed);

    % Define connectivity between RGC mosaic and cones
    [RGC_rf_map, RGC_rf_weights, extras] = cone_to_RGC_connectivity(cone_locations, RGC_locations, 'seed', 5555);            
    

    temp_std = zeros(num_cycles,1);
    for param_index = 1:length(deltapL)
        for cycle = 1:num_cycles

            if deltapL(param_index) == 0
                % random cone mosaic
                rand('twister', 1111+cycle);
                temp_cone_ids = rand(sqrt(num_cones));
                l_cone_ids = find(temp_cone_ids <= l_proportion);
                m_cone_ids = find(temp_cone_ids >= l_proportion);
                cone_types(l_cone_ids) = 0;
                cone_types(m_cone_ids) = 1;
            else
                clump_seed = 1111+cycle;
                [cone_types, extras] = filter_to_clump(cone_locations, filter_sigma, pL, deltapL(param_index), 'seed', clump_seed);
            end


            %artifactual clumping
%             combined_cone_fraction = 0.05;
%             num_combined_cones = round(num_cones * combined_cone_fraction);
%             temp_cone_dists = ipdm(cone_locations,'Subset', 'smallestfew', 'limit', num_combined_cones, 'result', 'struct');
%             flips = rand(num_combined_cones, 1);
%             num_changed_cones = 0;
%             for cone = 1:num_combined_cones
%                 if cone_types(temp_cone_dists.rowindex(cone)) ~= cone_types(temp_cone_dists.columnindex(cone))
%                     num_changed_cones = num_changed_cones + 1;
%                     if flips(cone) < 0.59
%                         cone_types(temp_cone_dists.rowindex(cone)) = 0;
%                         cone_types(temp_cone_dists.columnindex(cone)) = 0;
%                   else
%                         cone_types(temp_cone_dists.rowindex(cone)) = 1;
%                         cone_types(temp_cone_dists.columnindex(cone)) = 1;
%                     end
%                 end
%             end

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

            
            connectivity_matrix = zeros(num_cones, num_roi_RGCs);
            for RGC = 1:num_roi_RGCs
                temp_rf = zeros(num_cones,1);
                temp_rf(roi_RGC_rf_map{RGC}) = roi_RGC_rf_weights{RGC};
                connectivity_matrix(:,RGC) = temp_rf;
            end
            
            % normalize the area of the RF to 1
            if normalization
                for RGC = 1:num_roi_RGCs
                    center_cone_indices = find(connectivity_matrix(:,RGC) > 0);
                    total_center_input = sum(connectivity_matrix(center_cone_indices,RGC));
                    connectivity_matrix(:,RGC) = connectivity_matrix(:,RGC) ./ total_center_input;
                end
            end
            
            
            
            
            clear new_cone_types
            for cone = 1:num_cones
                if cone_types(cone) == 0
                    new_cone_types(cone) = 'L';
                else
                    new_cone_types(cone) = 'M';
                end
            end

            temp_indices = compute_opponency_index(connectivity_matrix, new_cone_types);
            temp_std(cycle) = std(temp_indices);
        end
        mean_oi(converg_index,param_index) = mean(temp_std);
        error_oi(converg_index,param_index) = std(temp_std)./sqrt(num_cycles);
    end
end    




    