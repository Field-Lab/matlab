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
%path_and_name{10,1} = '/snle/lab/Experiments/Array/Analysis/2007-03-27-2/data014/data014/data014';
%path_and_name{10,2} = 'grapes';
path_and_name{10,1} = '/snle/lab/Experiments/Array/Analysis/2008-08-27-0/data001-s5661-s9260/data001-s5661-s9260';
path_and_name{10,2} = 'peach';


num_datasets = length(path_and_name(:,1));
verbose = 0;

mosaic_counter = 0;
clear mosaic_purities mosaic_Gaussian_purities

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


        % get new coms and fits for cells
        % new coms
        ss_params.thresh = 3.5;
        distance_threshold = 40;
        temp_cell_indices = get_cell_indices(datarun, {cell_type});
        for RGC = 1:length(temp_cell_indices);
            temp_sta = get_sta(datarun, datarun.cell_ids(temp_cell_indices(RGC)));
            new_datarun.stas.stas{RGC} = temp_sta;
            datarun.stas.rf_coms{temp_cell_indices(RGC)} = rf_com(temp_sta,'ss_params',ss_params, 'distance_thresh', distance_threshold, 'positive', true);
        end
        % new fits
        test_datarun = fit_cone_rfs(datarun, {cell_type},'verbose', true, 'show_fit_params', false, 'foa_profile', [], 'foa_2d', []);





        % extract connectivity
        [mosaic_weights, selection, extras] = select_cone_weights(datarun, {cell_type}, 'thresh', 0.05, 'radius', [0 6], 'polarity', 1,...
                                                    'contiguity', true, 'scale', 3.0, 'remove_cones', 'SU');                                        
        connectivity = mosaic_weights .* selection;
        new_datarun = extras.new_datarun;
        [num_cones, num_RGCs] = size(connectivity);
        
        % normalize cone weights
        for RGC = 1:num_RGCs
            connectivity(:,RGC) = connectivity(:,RGC) ./ sum(connectivity(:,RGC));
        end

        num_bins = 1000;
        % get variability on weights
        ordered_cone_weights = cell(num_bins,1);
        convergences = zeros(num_RGCs,1);
        for RGC = 1:num_RGCs
            temp_fit = new_datarun.cones.rf_fits{RGC};
            roi_RGC_locations(RGC,:) = temp_fit.center;
            temp_cone_indices = find(connectivity(:,RGC) > 0);
            temp_cone_weights = connectivity(temp_cone_indices, RGC);
            temp_distances = ipdm(new_datarun.cones.centers(temp_cone_indices,:), temp_fit.center);
            temp_num_cones = length(temp_cone_weights);

            [sorted_distances, sorted_indices] = sort(temp_distances, 'ascend');
            sorted_weights = temp_cone_weights(sorted_indices);

            % store convergences
            convergences(RGC) = temp_num_cones;
            
            for cone = 1:temp_num_cones
                temp_weights = ordered_cone_weights{cone};
                temp_weights = [temp_weights, sorted_weights(cone)];
                ordered_cone_weights{cone} = temp_weights;
            end
        end

        clear weight_error
        for cone = 1:num_bins;
            temp_weights = ordered_cone_weights{cone};
            if ~isempty(temp_weights)
                weight_error(cone) = std(temp_weights);
            else
                weight_error(cone) = 0;
            end
        end
        mean_weight_error = mean(weight_error);
        
        figure(2)
        clf 
        hold on
        for cone = 1:num_bins
            temp_weights = ordered_cone_weights{cone};
            if isempty(temp_weights)
                temp_weights = 0;
            end
            errorbar(cone, mean(temp_weights), std(temp_weights), 'ko')
        end                                
                                        


        Gaussian_connectivity = zeros(num_cones, num_RGCs);
        for RGC = 1:num_RGCs
            temp_cone_indices = find(connectivity(:,RGC) > 0);
            temp_RF_center = new_datarun.cones.rf_fits{RGC}.center;
            temp_weights = connectivity(temp_cone_indices, RGC);


            %temp_distances = ipdm(temp_RF_center, new_datarun.cones.centers(temp_cone_indices,:));
            %temp_coef = [1, 5];
            %temp_fit_coef = nlinfit(temp_distances, temp_weights', 'zeroed_gaussian', temp_coef);
            %temp_Gaussian_weights = zeroed_gaussian(temp_fit_coef, temp_distances);

            temp_rf_fit = test_datarun.cones.rf_fit{temp_cell_indices(RGC)};
            

            temp_Gaussian_weights = 
 
%             figure(1)
%             clf
%             hold on
%             plot(temp_distances, temp_weights, 'ko')
%             plot(temp_distances, temp_Gaussian_weights, 'ro')
%             hold off
%             pause
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % get new center of mass for each RGC for fit

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            Gaussian_connectivity(temp_cone_indices, RGC) = temp_Gaussian_weights;
        end
        
        num_iters = 100;
        Gaussian_purity_sd = zeros(num_iters, 1);
        permuted_purity_sd = zeros(num_iters, 1);
        for iter = 1:num_iters
            new_Gaussian_connectivity = zeros(num_cones, num_RGCs);
            for RGC = 1:num_RGCs
                temp_cone_indices = find(Gaussian_connectivity(:, RGC) > 0);
                temp_num_cones = length(temp_cone_indices);
                temp_cone_weights = Gaussian_connectivity(temp_cone_indices,RGC);
                temp_weight_noise = 1.0*normrnd(0, weight_error(1:temp_num_cones));
                temp_cone_weights = abs(temp_cone_weights + temp_weight_noise');
                new_Gaussian_connectivity(temp_cone_indices, RGC) = temp_cone_weights;
            end
            
            rand_cone_indices = randperm(length(new_datarun.cones.types));
            permuted_purity_indices = compute_opponency_index(new_Gaussian_connectivity, new_datarun.cones.types(rand_cone_indices));
            permuted_purity_sd(iter) = std(permuted_purity_indices);
            
            Gaussian_purity_indices = compute_opponency_index(new_Gaussian_connectivity, new_datarun.cones.types);
            Gaussian_purity_sd(iter) = std(Gaussian_purity_indices);

            permuted_cone_PIs = compute_opponency_index(connectivity, new_datarun.cones.types(rand_cone_indices));
            permuted_cone_PI_SD(iter) = std(permuted_cone_PIs);
        end
        
        mosaic_counter = mosaic_counter +1;
        Gaussian_purity(mosaic_counter) = mean(Gaussian_purity_sd);
        Gaussian_purity_error(mosaic_counter) = std(Gaussian_purity_sd);
                
        permuted_Gaussian_purity(mosaic_counter) = mean(permuted_purity_sd);
        permuted_Gaussian_purity_error(mosaic_counter) = std(permuted_purity_sd);

        purity_indices = compute_opponency_index(connectivity, new_datarun.cones.types);
        purity_sd = std(purity_indices);
        mosaic_purities(mosaic_counter) = purity_sd
        
        permuted_cone_purity(mosaic_counter) = mean(permuted_cone_PI_SD);
        permuted_cone_purity_error(mosaic_counter) = std(permuted_cone_PI_SD);
     
        % permuted cones
%         rand_indices = randperm(length(connectivity_extras.cone_types));
%         permuted_purity_indices = compute_opponency_index(connectivity, connectivity_extras.cone_types(rand_indices));
%         permuted_purity_sd = std(permuted_purity_indices)
%         [permuted_purity_hist, hist_bins] = hist(permuted_purity_indices, [-1:0.1:1.0]);
% 
%         permuted_Gaussian_purity_indices = compute_opponency_index(Gaussian_connectivity, connectivity_extras.cone_types(rand_indices));
%         permuted_Gaussian_purity_sd = std(permuted_Gaussian_purity_indices)
%         [permuted_Gaussian_purity_hist, hist_bins] = hist(permuted_Gaussian_purity_indices, [-1:0.1:1.0]);


        
    end
end

figure(1)
clf
hold on
errorbar(mosaic_purities, Gaussian_purity, Gaussian_purity_error,'ko')
hold on
plot([0 0.5], [0 0.5], 'k')
axis([0 0.5 0 0.5])
axis square
xlabel('data')
ylabel('Gaussian w/ noise')
title('purity')
hold off
print(1, '/snle/home/gfield/Desktop/data-Gaussian','-dpdf')

figure(2)
clf
hold on
errorbar(permuted_cone_purity, permuted_Gaussian_purity, permuted_Gaussian_purity_error,'ko')
hold on
plot([0 0.5], [0 0.5], 'k')
axis([0 0.5 0 0.5])
axis square
xlabel('permuted')
ylabel('Gaussian w/ permuted')
title('purity')
hold off
print(2, '/snle/home/gfield/Desktop/permuted-Gaussian','-dpdf')


