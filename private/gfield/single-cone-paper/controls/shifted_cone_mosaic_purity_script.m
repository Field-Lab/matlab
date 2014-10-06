% path_and_name{1,1} = '/snle/lab/Experiments/Array/Analysis/2008-08-27-5/data003/data003/data003';
% path_and_name{1,2} = 'plantain';
% path_and_name{2,1} = '/snle/lab/Experiments/Array/Analysis/2008-08-26-2/data001-s6369-s9551/data001-s6369-s9551';
% path_and_name{2,2} = 'blueberry';
% path_and_name{3,1} = '/snle/lab/Experiments/Array/Analysis/2008-05-13-3/data006/data006';
% path_and_name{3,2} = 'kiwi';
% path_and_name{4,1} = '/snle/lab/Experiments/Array/Analysis/2008-08-27-0/data001-s5661-s9260/data001-s5661-s9260';
% path_and_name{4,2} = 'peach';
% 

% initialize counters
mosaic_counter = 0;
rf_center_threshold = 0.1;

% get paths to data of interest and number of datasets
[LMS_paths, LMS_names, cone_paths] = get_LMS_paths('high', 'cone_finding', 'standard');
num_datasets = length(LMS_paths);

for dataset = 1:num_datasets

    % load information from a data set
    clear datarun new_datarun sim_datarun cell_types
    datarun = load_data(LMS_paths{dataset});
    datarun = load_params(datarun,struct('verbose',1));  
    datarun = load_sta(datarun,'load_sta',[]);
    datarun = import_single_cone_data(datarun, LMS_names{dataset});    
    datarun.cones.mosaic = make_mosaic_struct(datarun.cones.centers);

    % get information on number of midget cells
    num_on_midgets = length(datarun.cell_types{3}.cell_ids);
    num_off_midgets = length(datarun.cell_types{4}.cell_ids);

    if num_on_midgets > 50
        on_midget_flag = true;
    else
        on_midget_flag = false;
    end
    if num_off_midgets > 50
        off_midget_flag = true;
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
    
    for tp = 1:length(cell_types)
        cell_type = cell_types(tp);
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get original connectivity matrix
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % extract connectivity
        [mosaic_weights, selection, extras] = select_cone_weights(datarun, {cell_type},...
                                                    'thresh', rf_center_threshold,...
                                                    'radius', [0 inf], 'polarity', 1,...
                                                    'contiguity', true,'scale', 3.0,...
                                                    'remove_cones', 'S');   
        

        connectivity = mosaic_weights .* selection; % keep weights continuous valued
        
        % the new datarun has S cone excluded from all fields of datarun.cone 
        new_datarun = extras.new_datarun;

        [num_cones, num_RGCs] = size(connectivity);
        %cell_num_tracker(mosaic_counter) = num_RGCs;



%         radius_scaler = 2.0;  % HARDCODED
%         RGC_min_convergence = 3;  % HARDCODED
%         RGC_max_convergence = 50;  % HARDCODED
% 
%         [connectivity, new_datarun, connectivity_extras] = extract_connectivity(datarun, cell_type, 'remove_cones', 'SU', 'required_sign', 'positive',...
%                                                     'min_convergence', RGC_min_convergence, 'max_convergence', RGC_max_convergence,...
%                                                     'max_radius', radius_scaler, 'normalize', true);
% 
%         num_roi_RGCs = length(connectivity(1,:));

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % shift cone mosaic
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get mean RF size
        RGC_indices = get_cell_indices(datarun, {cell_type});
        num_RGCs = length(RGC_indices);
        num_cones = length(datarun.cones.types);

        radii = zeros(num_RGCs,1);
        sorted_weights = cell(num_RGCs,1);

        for RGC = 1:num_RGCs
            temp_fit = new_datarun.cones.rf_fits{RGC};
            radii(RGC) = temp_fit.center_radius;
        end

        cone_shift_radius = median(radii) * 1.0;


        num_iter = 100;
        num_total_cones = length(new_datarun.cones.types);
        for iter = 1:num_iter
            rand_angle = 2*pi*rand(1);
            x_offset = cone_shift_radius * cos(rand_angle);
            y_offset = cone_shift_radius * sin(rand_angle);

            sim_datarun = new_datarun;
            sim_datarun.cones.centers = sim_datarun.cones.centers + repmat([x_offset, y_offset],num_total_cones, 1);

            sim_connectivity = zeros(num_total_cones, num_RGCs);
            for RGC = 1:num_RGCs;
                % get fit to use RGC center point
                temp_fit = new_datarun.cones.rf_fits{RGC};
                % get weights for RGC
                temp_cone_indices = find(connectivity(:,RGC)>0);
                % how many cones the RGC originally contacted
                temp_num_cones = length(temp_cone_indices);
                % cone weights of the RGC
                temp_weights = connectivity(temp_cone_indices, RGC);

                % order the distances from the RGC according to the cone distances
                % get distances to cones
                temp_cone_distances = ipdm(new_datarun.cones.centers(temp_cone_indices,:), temp_fit.center);
                % sort distances
                [sorted_distances, sorted_distance_indices] = sort(temp_cone_distances, 'ascend');
                % sort weights by distances
                distance_sorted_weights = temp_weights(sorted_distance_indices);

                % get the distance from RGC center point to the shifted cones
                temp_distances = ipdm(sim_datarun.cones.centers, temp_fit.center);
                % sort distances
                [sorted_distances, sorted_cone_indices] = sort(temp_distances, 'ascend');
                % assign the closest cones the RGC RF weights
                sim_connectivity(sorted_cone_indices(1:temp_num_cones), RGC) = distance_sorted_weights;
            end

            temp_indices = compute_opponency_index(sim_connectivity, sim_datarun.cones.types);
            shifted_cone_PIs(iter) = std(temp_indices);

        end

        mean_shifted_cones = mean(shifted_cone_PIs);
        sd_shifted_cones = std(shifted_cone_PIs);

 
        mosaic_counter = mosaic_counter + 1;
        shifted_cone_purity(mosaic_counter) = mean_shifted_cones;
        shifted_cone_purity_error(mosaic_counter) = sd_shifted_cones;

        purity_indices = compute_opponency_index(connectivity, new_datarun.cones.types);
        data_purity(mosaic_counter) = std(purity_indices);

    end
    
end

figure(1)
clf
hold on
errorbar(data_purity, shifted_cone_purity, shifted_cone_purity_error, 'ko')
plot([0 0.5], [0 0.5], 'k') 
axis([0 0.5 0 0.5])
axis square
xlabel('original cone locations')
ylabel('shifted cone locations')
title('purity index')
hold off
print(, '/snle/home/gfield/Desktop/fine-test-two','-dpdf')   



