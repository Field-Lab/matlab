clear 

% set analysis parameters
binarize_weights_flag = false; % if true, weights are binarized
use_cone_files = true; % if true, cone files are loaded from server
clumped_flag = false; % if true, clumped cone mosaics are loaded, 'use_cone_files', must also be true
rf_center_threshold = 0.1; % defines threshold for what is RF center
cell_cutoff = 50; % number of rgcs that are required to analyze the mosaic
num_permuted_cone_mosaics = 16;  % num
verbose = 0;


mosaic_counter = 0;

% get paths to data of interest and number of datasets
[LMS_paths, LMS_names, cone_paths] = get_LMS_paths('high', 'cone_finding', 'standard');
num_datasets = length(LMS_paths);

for dataset = 1:num_datasets

    % load information from a data set
    clear datarun new_datarun sim_datarun cell_types
    datarun = load_data(LMS_paths{dataset});
    datarun = load_params(datarun,struct('verbose',1));  
    datarun = load_sta(datarun,'load_sta',[]);
    datarun = import_single_cone_data(datarun, cone_paths{dataset});    
    datarun.cones.mosaic = make_mosaic_struct(datarun.cones.centers);

    % get information on number of midget cells
    num_on_midgets = length(datarun.cell_types{3}.cell_ids);
    num_off_midgets = length(datarun.cell_types{4}.cell_ids);


    % determine whether to analyze ON, OFF, or both midget mosaics
    if num_on_midgets > cell_cutoff
        on_midget_flag = true;
        %total_rgcs = total_rgcs + num_on_midgets;
    else
        on_midget_flag = false;
    end
    if num_off_midgets > cell_cutoff
        off_midget_flag = true;
        %total_rgcs = total_rgcs + num_off_midgets;
    else
        off_midget_flag = false;
    end
        
    if on_midget_flag && off_midget_flag
        cell_types = [3,4];
    elseif on_midget_flag && ~off_midget_flag
        cell_types = 3;
    elseif ~on_midget_flag && off_midget_flag
        cell_types = 4;
    else
        cell_types = [];
    end

  
    % analyze mosaic
    for tp = 1:length(cell_types)
        mosaic_counter = mosaic_counter + 1;
        cell_type = cell_types(tp);

        % get Gaussian weights for RGCs
        cell_indices = get_cell_indices(datarun, {cell_type});
        num_RGCs = length(cell_indices);
        
%% substitute Gaussian weights;
        for rgc = 1:num_RGCs
            temp_fit = datarun.cones.rf_fits{cell_indices(rgc)};
            temp_cone_locations = datarun.cones.centers;
            rgc_com = temp_fit.center;
            cone_distances = ipdm(rgc_com, temp_cone_locations);
            center_weights = temp_fit.center_scale * exp(-(cone_distances.^2) ./ (2*(temp_fit.center_radius).^2));
            surround_weights = temp_fit.surround_scale * exp(-(cone_distances.^2) ./ (2*(temp_fit.surround_radius).^2));
            %center_weights = temp_fit.center_scale * normpdf(cone_distances, 0, temp_fit.center_radius);
            %surround_weights = temp_fit.surround_scale * normpdf(cone_distances, 0, temp_fit.surround_radius);
            fit_weights = center_weights - surround_weights;
            datarun.cones.weights(:,cell_indices(rgc)) = fit_weights;
        end
        
%%
        
        % extract connectivity
        [mosaic_weights, selection, extras] = select_cone_weights(datarun, {cell_type},...
                                                    'thresh', rf_center_threshold,...
                                                    'radius', [0 6], 'polarity', 1,...
                                                    'contiguity', false,'scale', 3.0,...
                                                    'remove_cones', 'SU');   
        

        connectivity = mosaic_weights .* selection; % keep weights continuous valued
        
        %tmp_inds = find(connectivity(:,1));
        %tmp_wghts = connectivity(tmp_inds,1)
        
        % the new datarun has S cone excluded from all fields of datarun.cone 
        new_datarun = extras.new_datarun;

        [num_cones, num_RGCs] = size(connectivity);
        cell_num_tracker(mosaic_counter) = num_RGCs;

%%
%         % get Gaussian weights for RGCs
%         cell_indices = get_cell_indices(datarun, {cell_type});
%         new_weights = zeros(num_cones, num_RGCs);
%         for rgc = 1:num_RGCs
%             temp_fit = datarun.cones.rf_fits{cell_indices(rgc)};
%             temp_cone_indices = find(connectivity(:,rgc));
%             temp_cone_locations = datarun.cones.centers(temp_cone_indices,:);
%             rgc_com = temp_fit.center;
%             cone_distances = ipdm(rgc_com, temp_cone_locations);
%             center_weights = temp_fit.center_scale * normpdf(cone_distances, 0, temp_fit.center_radius);
%             surround_weights = temp_fit.surround_scale * normpdf(cone_distances, 0, temp_fit.surround_scale);
%             fit_weights = center_weights - surround_weights;
%             
%             if length(fit_weights) > 20;
%                 [sorted_weights, sorted_weight_indices] = sort(fit_weights, 'descend');
%                 new_weights(temp_cone_indices(sorted_weight_indices(1:20)), rgc) = fit_weights(sorted_weight_indices(1:20));
%             else
%                 new_weights(temp_cone_indices, rgc) = fit_weights;
%             end
%         end
%         
%         connectivity = new_weights;
%             
%%        
        % get pointers to sampled cones
        % purpose: this ensures that on cones sampled by the mosaic are permuted
        % this ensures that L:M cone ratio of the sampled cones is the same per and post permutation
        summed_selection = sum(selection, 2);
        cone_indices = find(summed_selection > 0);
        clear cone_pointers
        cone_pointers(1:num_cones) = false;
        cone_pointers(cone_indices) = true;

        % normalize cone weights
%         for RGC = 1:num_RGCs
%             connectivity(:,RGC) = connectivity(:,RGC) ./ sum(connectivity(:,RGC));
%         end

        % compute purity of cone mosaic
        purity_indices = compute_opponency_index(connectivity, new_datarun.cones.types);
        finited_indices = isfinite(purity_indices);
        purity_sd = std(purity_indices(finited_indices));
        mosaic_purities(mosaic_counter) = purity_sd;

        % determines wheter pre-calculated permutations are loaded from the server,
        % or whether permutations are calculated online.
        if use_cone_files  
            % load cone files
            for cone_mosaic = 1:num_permuted_cone_mosaics
                clear cleaned_cone_types temp_types cone_types
                suffix_string = num2str(11110+cone_mosaic);
                if clumped_flag
                    str2load = ['/snle/lab/Experiments/Array/Shared/one/simulations/clumped/',LMS_names{dataset},'-',suffix_string,'.txt'];
                else
                    str2load = ['/snle/lab/Experiments/Array/Shared/one/simulations/permuted/',LMS_names{dataset},'-',suffix_string,'.txt'];  
                end
                
                cones_file = dlmread([str2load],'\t');
                cone_types = cones_file(:,4);
                S_indices = find(cone_types == 3); % s cones
                U_indices = find(cone_types == 0); % untyped cones
                keep_cone_indices = setdiff([1:length(cone_types)], [S_indices; U_indices]);    
            
                temp_types = cone_types(keep_cone_indices);
                cleaned_cone_types(temp_types == 1,1) = 'L';
                cleaned_cone_types(temp_types == 2,1) = 'M';

                permuted_cone_PIs = compute_opponency_index(connectivity, cleaned_cone_types);
                finite_indices = isfinite(permuted_cone_PIs);
                permuted_cone_PI_SD(cone_mosaic) = std(permuted_cone_PIs(finite_indices));
     
                
                
                
%                 cones_file = dlmread([str2load],'\t');
%                 cone_types = cones_file(:,4);
%             
%                 cleaned_cone_types(cone_types == 1,1) = 'L';
%                 cleaned_cone_types(cone_types == 2,1) = 'M';
%                 cleaned_cone_types(cone_types == 3,1) = 'S';
%                 cleaned_cone_types(cone_types == 3,1) = 'U';
% 
%                 permuted_cone_PIs = compute_opponency_index(connectivity, cleaned_cone_types);
%                 finite_indices = isfinite(permuted_cone_PIs);
%                 permuted_cone_PI_SD(cone_mosaic) = std(permuted_cone_PIs(finite_indices));
            end
        else  
            % don't load cone files
            % permute cone labels within code -- clumping not supported
            for iter = 1:num_permuted_cone_mosaics
                new_cones = new_datarun.cones.types(cone_pointers);
                new_cone_num = length(new_cones);
                rand_cone_indices = randperm(new_cone_num);
                permuted_cone_PIs = compute_opponency_index(connectivity(cone_pointers,:), new_cones(rand_cone_indices));
 
                finite_indices = isfinite(permuted_cone_PIs);
                permuted_cone_PI_SD(iter) = std(permuted_cone_PIs(finite_indices));
            end
        end
        permuted_cone_purity(mosaic_counter) = mean(permuted_cone_PI_SD);
        permuted_cone_purity_error(mosaic_counter) = std(permuted_cone_PI_SD);

    end
end


figure
clf
hold on
errorbar(mosaic_purities, permuted_cone_purity, permuted_cone_purity_error,'ko')
hold on
plot([0 0.5], [0 0.5], 'k')
axis([0 0.5 0 0.5])
axis square
xlabel('data')
ylabel('permute')
title('purity')
hold off
print(6, '/snle/home/gfield/Desktop/gauss_weights.pdf','-dpdf')



