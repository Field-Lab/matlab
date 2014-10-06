path_and_name{1,1} = '/snle/lab/Experiments/Array/Analysis/2008-08-27-5/data003/data003/data003';
path_and_name{1,2} = 'plantain';
path_and_name{2,1} = '/snle/lab/Experiments/Array/Analysis/2008-08-26-2/data001-s6369-s9551/data001-s6369-s9551';
path_and_name{2,2} = 'blueberry';
path_and_name{3,1} = '/snle/lab/Experiments/Array/Analysis/2008-05-13-3/data006/data006';
path_and_name{3,2} = 'kiwi';
path_and_name{4,1} = '/snle/lab/Experiments/Array/Analysis/2008-08-27-0/data001-s5661-s9260/data001-s5661-s9260';
path_and_name{4,2} = 'peach';


num_datasets = length(path_and_name(:,1));
verbose = 0;

mosaic_counter = 0;
clear Gaussian_purity shifted_cone_purity shifted_cone_purity_error
%weight_method = 'Gaussian'; % or 'shuffle'
weight_method = 'shuffle'; % or 'shuffle'

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

        num_roi_RGCs = length(connectivity(1,:));
        ordered_cone_weights = cell(RGC_max_convergence,1);
        if strcmp(weight_method,'shuffle')
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
        end
                                                

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % shift cone mosaic
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get mean RF size
        RGC_indices = get_cell_indices(datarun, {cell_type});
        num_RGCs = length(RGC_indices);
        num_cones = length(datarun.cones.types);

        radii = zeros(num_roi_RGCs,1);
        sorted_weights = cell(num_roi_RGCs,1);
        clear fit_coefs
        for RGC = 1:num_roi_RGCs
            temp_fit = new_datarun.cones.rf_fits{RGC};
            radii(RGC) = temp_fit.center_radius;
            scales(RGC) = temp_fit.center_scale;            
            
            if strcmp(weight_method, 'Gaussian')
                temp_cone_indices = find(connectivity(:,RGC) > 0);
                temp_weights = connectivity(temp_cone_indices, RGC);

                temp_distances = ipdm(new_datarun.cones.centers(temp_cone_indices,:), temp_fit.center);
                temp_coef = [0.2, 5];
                temp_fit_coef = nlinfit(temp_distances, temp_weights, 'zeroed_gaussian', temp_coef);
                temp_cone_weights = zeroed_gaussian(temp_fit_coef, temp_distances);

                [temp_sorted_dists temp_sorted_indices] = sort(temp_distances, 'ascend');
                temp_cone_weights = temp_cone_weights(temp_sorted_indices);
                sorted_weights{RGC} = temp_cone_weights;
            end
            if strcmp(weight_method, 'shuffle')
                temp_cone_indices = find(connectivity(:,RGC) > 0);
                temp_num_cones = length(temp_cone_indices);
                new_weights = zeros(temp_num_cones,1);
                for cone = 1:temp_num_cones
                    temp_weights = ordered_cone_weights{cone};
                    rand_index = floor(rand(1) * length(temp_weights)) + 1;
                    new_weights(cone) = temp_weights(rand_index);
                end
                sorted_weights{RGC} = new_weights;
            end
            if strcmp(weight_method, 'static')
                temp_cone_indices = find(connectivity(:,RGC) > 0);
                temp_num_cones = length(temp_Cone_indices);
                new_weights = 
                
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

            sim_connectivity = zeros(num_total_cones, num_roi_RGCs);
            for RGC = 1:num_roi_RGCs;
                temp_weights = sorted_weights{RGC};
                temp_fit = sim_datarun.cones.rf_fits{RGC};
                temp_num_cones = length(temp_weights);

                temp_distances = ipdm(sim_datarun.cones.centers, temp_fit.center);
                [sorted_distances, sorted_cone_indices] = sort(temp_distances, 'ascend');
                sim_connectivity(sorted_cone_indices(1:temp_num_cones), RGC) = temp_weights;
            end

            temp_indices = compute_opponency_index(sim_connectivity, sim_datarun.cones.types);
            shifted_cone_PIs(iter) = std(temp_indices);

        end

        mean_shifted_cones = mean(shifted_cone_PIs);
        sd_shifted_cones = std(shifted_cone_PIs);


        [shifted_cone_hist, hist_bins] = hist(shifted_cone_PIs, [0.2:0.01:0.6]);
        figure(1)
        bar(hist_bins, shifted_cone_hist, 'm')

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Gaussian weights
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        [num_cones, num_RGCs] = size(connectivity);                                      
        LM_cone_locations = new_datarun.cones.centers;                                         
        new_connectivity = zeros(num_cones, num_RGCs);

        for RGC = 1:num_RGCs
            temp_cone_indices = find(connectivity(:,RGC) > 0);
            temp_fit = new_datarun.cones.rf_fits{RGC};
            temp_RF_center = temp_fit.center;
            temp_distances = ipdm(temp_RF_center, LM_cone_locations(temp_cone_indices,:));
            if strcmp(weight_method, 'Gaussian')

                temp_coef = [0.2, 5];
                temp_fit_coef = nlinfit(temp_distances, connectivity(temp_cone_indices, RGC)', 'zeroed_gaussian', temp_coef);

                temp_cone_weights = zeroed_gaussian(temp_fit_coef, temp_distances);
                new_connectivity(temp_cone_indices, RGC) = temp_cone_weights;

            end
            if strcmp(weight_method, 'shuffle')
                [sorted_distances, sorted_indices] = sort(temp_distances, 'ascend');
                temp_weights = sorted_weights{RGC};
                new_connectivity(temp_cone_indices(sorted_indices), RGC) = temp_weights;
            end
        end

        new_purity_indices = compute_opponency_index(new_connectivity, new_datarun.cones.types);
        new_purity_sd = std(new_purity_indices);
        [new_purity_hist, hist_bins] = hist(new_purity_sd, [-1:0.1:1.0]);

        mosaic_counter = mosaic_counter + 1;
        new_purity(mosaic_counter) = new_purity_sd;
        shifted_cone_purity(mosaic_counter) = mean_shifted_cones;
        shifted_cone_purity_error(mosaic_counter) = sd_shifted_cones;

        purity_indices = compute_opponency_index(connectivity, new_datarun.cones.types);
        data_purity(mosaic_counter) = std(purity_indices);

    end
    
end

figure(1)
clf
hold on
errorbar(new_purity, shifted_cone_purity, shifted_cone_purity_error, 'ko')
plot([0 0.5], [0 0.5], 'k') 
axis([0 0.5 0 0.5])
xlabel('original cone locations')
ylabel('shifted cone locations')
title('purity: permuted weights')
hold off
print(5, '/snle/home/gfield/Desktop/fine-test-two','-dpdf')   

figure(2)
clf
hold on
plot(data_purity, new_purity, 'ko')
plot([0 0.5], [0 0.5], 'k') 
axis([0 0.5 0 0.5])
xlabel('data')
ylabel('shuffled')
title('purity')
hold off


