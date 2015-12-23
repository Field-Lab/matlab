[LMS_paths, LMS_names] = get_LMS_paths('high');

num_datasets = length(LMS_paths);
verbose = 0;

mosaic_counter = 0;

for dataset = 1:num_datasets
    clear datarun new_datarun sim_datarun
    datarun = load_data(LMS_paths{dataset});
    datarun = load_params(datarun,struct('verbose',1));  
    datarun = load_sta(datarun,'load_sta',[]);
    datarun = import_single_cone_data(datarun, LMS_names{dataset});
    
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
        mosaic_counter = mosaic_counter + 1;
        cell_type = cell_types(tp);
        
        temp_filename = [LMS_names{dataset},'-',num2str(tp)];
        load(temp_filename)

        temp_cell_indices = get_cell_indices(datarun, {cell_type});

        cell_numbers(mosaic_counter) = length(temp_cell_indices);

        % extract connectivity
        [mosaic_weights, selection, extras] = select_cone_weights(datarun, {cell_type}, 'thresh', 0.05, 'radius', [0 inf], 'polarity', 1,...
                                                    'contiguity', true, 'scale', 3.0, 'remove_cones', 'SU');   
        connectivity = mosaic_weights .* selection;
        new_datarun = extras.new_datarun;
        [num_cones, num_RGCs] = size(connectivity);   

       % replace weights with circular Gaussian weights;
        gauss_connectivity = zeros(num_cones, num_RGCs);
        for rgc = 1:num_RGCs
            cone_indices = find(connectivity(:,rgc) > 0);
            cone_centers = new_datarun.cones.centers(cone_indices,:);
            if isempty(cone_indices)
                rgc = rgc
                warning('RGC has no significant cones')
            end
            rgc_fit = new_datarun.cones.rf_fits{rgc};
            rgc_center = rgc_fit.center;
            cone_distances = ipdm(rgc_center, cone_centers);

            new_weights = normpdf(cone_distances, 0, rgc_fit.center_radius);
            
            %figure(1)
            %clf; hold on
            %plot(cone_distances, new_weights, 'ko')
            %plot(cone_distances, connectivity(cone_indices, rgc)./sum(connectivity(cone_indices, rgc)), 'ro')
            %hold off

            gauss_connectivity(cone_indices, rgc) = new_weights;

        end
        
       % calculate the PI for the new gaussian weights
        purity_indices = compute_opponency_index(gauss_connectivity, new_datarun.cones.types);
        purity_sd = std(purity_indices);
        gauss_purities(mosaic_counter) = purity_sd;

        % shuffle the new Gaussian weights
        num_iter = 50;
        sd_shuffled_purity = zeros(1, num_iter);
        rand('twister', 11111)
        for iter = 1:num_iter
            shuffled_gauss_connectivity = zeros(num_cones, num_RGCs);
            for rgc = 1:num_RGCs
                cone_indices = find(gauss_connectivity(:,rgc) > 0);
                temp_weights = gauss_connectivity(cone_indices, rgc);
                shuffled_indices = randperm(length(cone_indices));
                shuffled_weights = temp_weights(shuffled_indices);
                
                shuffled_gauss_connectivity(cone_indices,rgc) = shuffled_weights;
            end

            temp_purity_indices = compute_opponency_index(shuffled_gauss_connectivity, new_datarun.cones.types);
            sd_shuffled_purity(iter) = std(temp_purity_indices);
        end
        
        shuffled_gauss_purities(mosaic_counter) = mean(sd_shuffled_purity);
        shuffled_gauss_purities_error(mosaic_counter) = std(sd_shuffled_purity);

    end
end

figure(1)
clf
hold on
errorbar(gauss_purities, shuffled_gauss_purities, shuffled_gauss_purities_error, 'ko')
plot([0 0.5], [0 0.5], 'k-')



