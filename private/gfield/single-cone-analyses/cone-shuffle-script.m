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
path_and_name{10,1} = '/snle/lab/Experiments/Array/Analysis/2008-08-27-0/data001-s5661-s9260/data001-s5661-s9260';
path_and_name{10,2} = 'peach';


num_datasets = length(path_and_name(:,1));
verbose = 0;

mosaic_counter = 0;
clear mosaic_unique_cone_fraction mosaic_unique_cone_fraction_sd data_SD_PIs

rand('twister', 11111);
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

        
        %[mosaic_weights, selection] = select_cone_weights(datarun, {cell_type});                                        
                                                
        %connectivity = mosaic_weights .* selection;
        
        [num_cones, num_RGCs] = size(connectivity);
        
        num_shuffles = 100;
        shuffled_cones_SD_PIs = zeros(1,num_shuffles);
        for shuffle = 1:num_shuffles
            temp_datarun = new_datarun;
            fraction_shared = zeros(num_RGCs, 1);
            unique_cone_fraction = zeros(num_RGCs, 1);
            for RGC = 1:num_RGCs
                % find cones connectied to RGC of interest
                temp_cone_indices = find(connectivity(:,RGC) > 0);
                temp_num_center_cones = length(temp_cone_indices);
                % identify whether any of these cones are connected to another
                % RGC
                clear shared_cone_indices
                shared_cone_counter = 0;
                shared_cone_indices = [];
                for cone = 1:temp_num_center_cones
                    temp_test = find(connectivity(temp_cone_indices(cone), :) > 0);
                    if length(temp_test) > 5
                        shared_cone_counter = shared_cone_counter +1;
                        shared_cone_indices(shared_cone_counter) = temp_cone_indices(cone);
                    end

                end
                if isempty(shared_cone_indices);
                    exclusive_cones = temp_cone_indices;
                    fraction_shared(RGC) = 0;
                else
                    exclusive_cones = setdiff(temp_cone_indices, shared_cone_indices);
                    fraction_shared(RGC) = length(shared_cone_indices) ./ length(temp_cone_indices);
                end
  
                temp_RF_fits = new_datarun.cones.rf_fits{RGC};
                
                % visualize cone locations and mark shared cones
%                 figure(1)
%                 clf
%                 hold on
%                 plot(new_datarun.cones.centers(shared_cone_indices, 1), new_datarun.cones.centers(shared_cone_indices, 2), 'oc', 'MarkerFaceColor', 'c')
%                 plot(new_datarun.cones.centers(exclusive_cones, 1), new_datarun.cones.centers(exclusive_cones, 2), 'ok', 'MarkerFaceColor', 'k')
%                 plot(temp_RF_fits.center(1), temp_RF_fits.center(2), 'mo', 'MarkerFaceColor', 'm')
%                 hold off
%                 pause
                
                % permute cone labels accordingly
                num_excl_cones = length(exclusive_cones);
                rand_indices = randperm(num_excl_cones);
                temp_datarun.cones.types(exclusive_cones) = new_datarun.cones.types(exclusive_cones(rand_indices));
                
                % extract the fraction of RF center volume that is occupied
                % by unique cones
                total_center_volume = sum(connectivity(temp_cone_indices, RGC));
                unique_center_volume = sum(connectivity(exclusive_cones, RGC));
                unique_cone_fraction(RGC) = unique_center_volume ./ total_center_volume;
                
            end

            %simulations
            shuffled_cones_PIs = compute_opponency_index(connectivity, temp_datarun.cones.types);
            shuffled_cones_SD_PIs(shuffle) = std(shuffled_cones_PIs);
        end
        
        mosaic_counter = mosaic_counter + 1;

        %data
        data_PIs = compute_opponency_index(connectivity, new_datarun.cones.types);
        data_SD_PIs(mosaic_counter) = std(data_PIs);
        track_fraction_shared(mosaic_counter) = mean(fraction_shared);
        track_fraction_shared_sd(mosaic_counter) = std(fraction_shared);
        track_fraction_shared_error(mosaic_counter) = std(fraction_shared) ./ (sqrt(num_RGCs) - 1);
        
        mean_shuffled_cones_PI(mosaic_counter) = mean(shuffled_cones_SD_PIs);
        sd_shuffled_cones_PI(mosaic_counter) = std(shuffled_cones_SD_PIs);
        
        % store the average and std of the unique cone fraction
        mosaic_unique_cone_fraction(mosaic_counter) = mean(unique_cone_fraction);
        mosaic_unique_cone_fraction_sem(mosaic_counter) = std(unique_cone_fraction) ./ sqrt(num_RGCs);

    end
    
end

figure(20)
clf
hold on
errorbar(data_SD_PIs, mean_shuffled_cones_PI, sd_shuffled_cones_PI, 'ko')    
plot([0 0.5], [0 0.5], 'k') 
axis([0.2 0.5 0.2 0.5])
xlabel('data')
ylabel('cones locally shuffled')
title('purity')
axis square
hold off
%print(10, '/snle/home/gfield/Desktop/locally_shifted','-dpdf')   

   

figure(11)
clf
errorbar(track_fraction_shared, track_fraction_shared_error, 'ok')
xlabel('preparation')
ylabel('fraction of cones shared')
title('on and off midgets')
axis([0 20 0 1])

figure(12)
clf
errorbar(mosaic_unique_cone_fraction, mosaic_unique_cone_fraction_sem, 'ko')
xlabel('preparation')
ylabel('unique RF fraction')
axis([0 20 0 1])

