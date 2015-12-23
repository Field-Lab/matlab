path_and_name{1,1} = '/snle/lab/Experiments/Array/Analysis/2008-08-27-5/data003/data003/data003';
path_and_name{1,2} = 'plantain';
path_and_name{2,1} = '/snle/lab/Experiments/Array/Analysis/2008-04-22-5/data006/data006';
path_and_name{2,2} = 'plum';
path_and_name{3,1} = '/snle/lab/Experiments/Array/Analysis/2008-08-26-2/data001-s6369-s9551/data001-s6369-s9551';
path_and_name{3,2} = 'blueberry';
path_and_name{4,1} = '/snle/lab/Experiments/Array/Analysis/2008-12-12-1/data005/data005';
path_and_name{4,2} = 'butterfly';
path_and_name{5,1} = '/snle/lab/Experiments/Array/Analysis/2007-08-21-1/data003/data003';
path_and_name{5,2} = 'pomegranate';
path_and_name{6,1} = '/snle/lab/Experiments/Array/Analysis/2009-02-28-0/data006/data006';
path_and_name{6,2} = 'cherry';
path_and_name{7,1} = '/snle/lab/Experiments/Array/Analysis/2008-03-25-3/data002/data002';
path_and_name{7,2} = 'cherimoya';
path_and_name{8,1} = '/snle/lab/Experiments/Array/Analysis/2008-05-13-3/data006/data006';
path_and_name{8,2} = 'kiwi';
path_and_name{9,1} = '/snle/lab/Experiments/Array/Analysis/2008-04-30-2/data004/data004/data004';
path_and_name{9,2} = 'mango';
path_and_name{10,1} = '/snle/lab/Experiments/Array/Analysis/2007-03-27-2/data014/data014/data014';
path_and_name{10,2} = 'grapes';
path_and_name{11,1} = '/snle/lab/Experiments/Array/Analysis/2008-08-27-0/data001-s5661-s9260/data001-s5661-s9260';
path_and_name{11,2} = 'peach';


num_datasets = length(path_and_name(:,1));
verbose = 0;

mosaic_counter = 0;
clear permute_weights_purities permute_weights_purities_error data_purities

for dataset = 1:num_datasets
    clear datarun new_datarun sim_datarun
    datarun = load_data(path_and_name{dataset,1});
    datarun = load_params(datarun,struct('verbose',1));  
    datarun = load_sta(datarun,'load_sta',[]);
    datarun = import_single_cone_data(datarun, path_and_name{dataset, 2});
    
    num_on_midgets = length(datarun.cell_types{3}.cell_ids);
    num_off_midgets = length(datarun.cell_types{4}.cell_ids);
    
    if num_on_midgets > 20
        on_midget_flag = true;
    else
        on_midget_flag = false;
    end
    if num_off_midgets > 20
        off_midget_flag = true;
    else
        off_midget_flag = false;
    end
    
    if on_midget_flag && off_midget_flag
        cell_types = [3,4];
    elseif on_midget_flag && ~off_midget_flag
        cell_types = 3;
    elseif ~on_midget_flag && off_midget_flag
        cell_type = 4;
    else
        cell_types = [];
    end

    for tp = 1:length(cell_types)
        cell_type = cell_types(tp);
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get original connectivity matrix
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        radius_scaler = 2.0;  % HARDCODED
        RGC_min_convergence = 3;  % HARDCODED
        RGC_max_convergence = 50;  % HARDCODED

        [connectivity, new_datarun, connectivity_extras] = extract_connectivity(datarun, cell_type, 'remove_cones', 'SU', 'required_sign', 'positive',...
                                                    'min_convergence', RGC_min_convergence, 'max_convergence', RGC_max_convergence,...
                                                    'max_radius', radius_scaler, 'normalize', true);


                                                
        [num_roi_cones, num_roi_RGCs] = size(connectivity);
        ordered_cone_weights = cell(RGC_max_convergence,1);
        for RGC = 1:num_roi_RGCs
            temp_fit = new_datarun.cones.rf_fits{RGC};
            roi_RGC_locations(RGC,:) = temp_fit.center;
            temp_cone_indices = find(connectivity(:,RGC) > 0);
            temp_cone_weights = connectivity(temp_cone_indices, RGC);
            temp_distances = ipdm(new_datarun.cones.centers(temp_cone_indices,:), temp_fit.center);
            temp_num_cones = length(temp_cone_weights);

            [sorted_distances, sorted_indices] = sort(temp_distances, 'ascend');
            sorted_weights = temp_cone_weights(sorted_indices);

            for cone = 1:temp_num_cones
                temp_weights = ordered_cone_weights{cone};
                temp_weights = [temp_weights, sorted_weights(cone)];
                ordered_cone_weights{cone} = temp_weights;
            end
        end                                               

%         figure(2)
%         clf 
%         hold on
%         for cone = 1:RGC_max_convergence
%             temp_weights = ordered_cone_weights{cone};
%             if isempty(temp_weights)
%                 temp_weights = 0;
%             end
%             errorbar(cone, mean(temp_weights), std(temp_weights), 'ko')
%         end
        
        
        % Get cell_array of ordered cone weights
        distance_matrix = ipdm(new_datarun.cones.centers, roi_RGC_locations);

        
        %%% across cell weight shuffle
        num_shuffles = 100;
        for shuffle = 1:num_shuffles
            shuffled_connectivity = zeros(num_roi_cones, num_roi_RGCs);
            for RGC = 1:num_roi_RGCs
                temp_cone_indices = find(connectivity(:,RGC) > 0);
                temp_num_cones = length(temp_cone_indices);
                temp_distances = distance_matrix(temp_cone_indices, RGC);
                [sorted_distances, sorted_indices] = sort(temp_distances, 'ascend');
                
                new_weights = zeros(temp_num_cones,1);
                for cone = 1:temp_num_cones
                    temp_weights = ordered_cone_weights{cone};
                    rand_index = floor(rand(1) * length(temp_weights)) + 1;
                    new_weights(sorted_indices(cone)) = temp_weights(rand_index);
                end
                shuffled_connectivity(temp_cone_indices, RGC) = new_weights ./ sum(new_weights);
            end
            % compute opponency indices
            temp_indices = compute_opponency_index(shuffled_connectivity, new_datarun.cones.types);
            temp_sd = std(temp_indices);
            weight_permutation_PI_SD(shuffle) = temp_sd;

            % permute cones 
            rand_cone_indices = randperm(num_roi_cones);
            % w/ permuted weights
            temp_indices = compute_opponency_index(shuffled_connectivity, new_datarun.cones.types(rand_cone_indices));
            weight_and_cone_perm_PI_SD(shuffle) = std(temp_indices);
            % w/ original weights
            temp_indices = compute_opponency_index(connectivity, new_datarun.cones.types(rand_cone_indices));
            cone_perm_PI_SD(shuffle)= std(temp_indices);
            
        end
        
        mosaic_counter = mosaic_counter + 1;
        permute_weights_purities(mosaic_counter) = mean(weight_permutation_PI_SD);
        permute_weights_purities_error(mosaic_counter) = std(weight_permutation_PI_SD);
        permute_weights_and_cones(mosaic_counter) = mean(weight_and_cone_perm_PI_SD);
        permute_weights_and_cones_error(mosaic_counter) = std(weight_and_cone_perm_PI_SD);
        permute_cones_purity(mosaic_counter) = mean(cone_perm_PI_SD);
        permute_cones_purity_error(mosaic_counter) = std(cone_perm_PI_SD);
        
        %%% within cell weight shuffle
        num_shuffles = 100;
        within_cell_shuffle_purity = zeros(1,num_shuffles);
        for shuffle = 1:num_shuffles;
            shuffled_connectivity_matrix = zeros(num_roi_cones, num_roi_RGCs);
            for RGC = 1:num_roi_RGCs
                % get the weights for the RGC
                temp_cone_indices = find(connectivity(:, RGC) > 0);
                temp_cone_weights = connectivity(temp_cone_indices, RGC);
                new_weight_order = randperm(length(temp_cone_indices));

                shuffled_connectivity_matrix(temp_cone_indices, RGC) = temp_cone_weights(new_weight_order);
            end

            temp_indices = compute_opponency_index(shuffled_connectivity_matrix, new_datarun.cones.types);
            current_std = std(temp_indices);
            within_cell_shuffle_purity(shuffle) = current_std;
        end

        purity_within_cell_shuffle(mosaic_counter) = mean(within_cell_shuffle_purity);
        purity_within_cell_shuffle_error(mosaic_counter) = std(within_cell_shuffle_purity);
        
        temp_purities = compute_opponency_index(connectivity, new_datarun.cones.types);
        data_purities(mosaic_counter) = std(temp_purities);
        
       
    end

end



figure(1)
clf
hold on
errorbar(data_purities, permute_weights_purities, permute_weights_purities_error, 'ko')
plot([0 0.5], [0 0.5], 'k')
xlabel('data')
ylabel('permuted weights (across cells)')
title('purity')
axis([0 0.5 0 0.5])
axis square
hold off
print(1, '/snle/home/gfield/Desktop/permute_weights','-dpdf')   

figure(2)
clf
hold on
errorbar(permute_cones_purity, permute_weights_and_cones, permute_weights_and_cones_error, 'ko')
plot([0 0.5], [0 0.5], 'k')
xlabel('permuted cones')
ylabel('permuted cones and weights')
title('purity')
axis([0 0.5 0 0.5])
axis square
hold off
print(2, '/snle/home/gfield/Desktop/permute_weights_and_cones','-dpdf')   


