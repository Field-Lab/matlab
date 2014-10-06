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
clear data_PIs_SD shuffled_weights_mean shuffled_weights_error perm_cones_purity_mean perm_cone_weights_purity_mean
clear perm_cone_weights_purity_error perm_cone_purity_error shuffled_weights_perm_cones shuffled_weights_perm_cones_error

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
        %RGC_min_convergence = 3;  % HARDCODED
        %RGC_max_convergence = 50;  % HARDCODED

        %[connectivity, new_datarun, connectivity_extras] = extract_connectivity(datarun, cell_type, 'remove_cones', 'SU', 'required_sign', 'positive',...
        %                                            'min_convergence', RGC_min_convergence, 'max_convergence', RGC_max_convergence,...
        %                                            'max_radius', radius_scaler, 'normalize', true);


        
        
        
        % get SD(PI) of data
        data_PIs = compute_opponency_index(connectivity, new_datarun.cones.types);                                        
        mosaic_counter = mosaic_counter + 1;
        data_PIs_SD(mosaic_counter) = std(data_PIs);
        
        % get SD(PI) of shuffled_weights
        [num_roi_cones, num_roi_RGCs] = size(connectivity);                            
        num_shuffles = 100;
        shuffled_weights_PI_SD = zeros(num_shuffles,1);
        perm_weights_cones_PI_SD = zeros(num_shuffles,1);
        for shuffle = 1:num_shuffles
            shuffled_connectivity = zeros(num_roi_cones, num_roi_RGCs);
            for RGC = 1:num_roi_RGCs
                temp_cone_indices = find(connectivity(:,RGC) ~= 0);
                temp_num_cones = length(temp_cone_indices);
                rand_indices = randperm(temp_num_cones);
                shuffled_connectivity(temp_cone_indices, RGC) = connectivity(temp_cone_indices(rand_indices), RGC);
            end
            shuffled_weights_PIs = compute_opponency_index(shuffled_connectivity, new_datarun.cones.types);
            shuffled_weights_PI_SD(shuffle) = std(shuffled_weights_PIs);
            
            % permuted the cones after permuting the weights for additional
            % effects
            rand_cone_indices = randperm(num_roi_cones);
            perm_weights_cones_PIs = compute_opponency_index(shuffled_connectivity, new_datarun.cones.types(rand_cone_indices));          
            perm_weights_cones_PI_SD(shuffle) = std(perm_weights_cones_PIs);
            
            perm_again_indices = randperm(num_roi_cones);
            perm_again_PIs = compute_opponency_index(shuffled_connectivity, new_datarun.cones.types(rand_cone_indices(perm_again_indices)));
            perm_again_PI_SD(shuffle) = std(perm_again_PIs);
        end

        shuffled_weights_mean(mosaic_counter) = mean(shuffled_weights_PI_SD);
        shuffled_weights_error(mosaic_counter) = std(shuffled_weights_PI_SD);
        
        shuffled_weights_perm_cones(mosaic_counter) = mean(perm_weights_cones_PI_SD);
        shuffled_weights_perm_cones_error(mosaic_counter) = std(perm_weights_cones_PI_SD);
        
        perm_again_purity(mosaic_counter) = mean(perm_again_PI_SD);
        perm_again_purity_error(mosaic_counter) = std(perm_again_PI_SD);
        
        
        % get SD(PI) for permuted cones and weights
        num_shuffles = 100;
        perm_cones_PI_SD = zeros(num_shuffles, 1);
        perm_cones_weights_PI_SD = zeros(num_shuffles, 1);

        
        for shuffle = 1:num_shuffles
            perm_cones_weights_connectivity = zeros(num_roi_cones, num_roi_RGCs);
            
            rand_cone_indices = randperm(num_roi_cones);
            perm_cones_PIs = compute_opponency_index(connectivity, new_datarun.cones.types(rand_cone_indices));
            perm_cones_PI_SD(shuffle) = std(perm_cones_PIs);
            
            for RGC = 1:num_roi_RGCs
                temp_cone_indices = find(connectivity(:,RGC) ~= 0);
                temp_num_cones = length(temp_cone_indices);
                rand_indices = randperm(temp_num_cones);
                perm_cones_weights_connectivity(temp_cone_indices, RGC) = connectivity(temp_cone_indices(rand_indices), RGC);
            end    
            perm_cones_weights_PIs = compute_opponency_index(perm_cones_weights_connectivity, new_datarun.cones.types(rand_cone_indices));
            perm_cones_weights_PI_SD(shuffle) = std(perm_cones_weights_PIs);
    
        end
        
        perm_cones_purity_mean(mosaic_counter) = mean(perm_cones_PI_SD);
        perm_cone_purity_error(mosaic_counter) = std(perm_cones_PI_SD);
        perm_cone_weights_purity_mean(mosaic_counter) = mean(perm_cones_weights_PI_SD);
        perm_cone_weights_purity_error(mosaic_counter) = std(perm_cones_weights_PI_SD);
    end
end


% compute z-score effects
total_zscores = (data_PIs_SD - perm_cones_purity_mean) ./ perm_cone_purity_error;
fine_zscores = (data_PIs_SD - shuffled_weights_mean) ./ shuffled_weights_error;


figure(5)
clf
hold on
plot(total_zscores, fine_zscores, 'ko')
plot([0 10], [0 10], 'k') 
axis([0 10 0 10])
axis square
xlabel('total')
ylabel('fine')
hold off
print(5, '/snle/home/gfield/Desktop/zscores','-dpdf') 

figure(1)
clf
hold on
errorbar(data_PIs_SD, shuffled_weights_mean, shuffled_weights_error, 'ko')
plot([0 0.5], [0 0.5], 'k') 
axis([0 0.5 0 0.5])
axis square
xlabel('data')
ylabel('permuted weights (within cells)')
title('purity')
hold off
print(1, '/snle/home/gfield/Desktop/perm_weights','-dpdf') 

figure(2)
clf
hold on
errorbar(shuffled_weights_mean, shuffled_weights_perm_cones, shuffled_weights_error, 'ko')
plot([0 0.5], [0 0.5], 'k') 
axis([0 0.5 0 0.5])
axis square
xlabel('permuted weights')
ylabel('permuted weights and cones')
title('purity')
hold off
print(2, '/snle/home/gfield/Desktop/perm_weights_cones','-dpdf') 


figure(3)
clf
hold on
errorbar(perm_cones_purity_mean, perm_cone_weights_purity_mean, perm_cone_weights_purity_error, 'ko')
plot([0 0.5], [0 0.5], 'k') 
axis([0 0.5 0 0.5])
axis square
xlabel('permuted cones')
ylabel('permuted weights and cones')
title('purity')
hold off
print(3, '/snle/home/gfield/Desktop/check','-dpdf')             
            
            

figure(4)
clf
hold on
errorbar(data_PIs_SD, perm_cones_purity_mean, perm_cone_purity_error, 'ko')
plot([0 0.5], [0 0.5], 'k') 
axis([0 0.5 0 0.5])
axis square
xlabel('data')
ylabel('permuted cones')
title('purity')
hold off
print(4, '/snle/home/gfield/Desktop/perm-cones','-dpdf')             
            
            
figure(6)
clf
hold on
errorbar(data_PIs_SD, shuffled_weights_perm_cones, shuffled_weights_perm_cones_error, 'ko')
errorbar(data_PIs_SD, shuffled_weights_mean, shuffled_weights_error, 'k.', 'MarkerSize', 18)
plot([0 0.5], [0 0.5], 'k') 
axis([0.2 0.5 0.2 0.5])
axis square
xlabel('permuted cones')
ylabel('permuted weights and cones')
title('purity')
hold off
print(6, '/snle/home/gfield/Desktop/check','-dpdf')             
            
            
            
        
        
        
        
        
                                            
                                                
                                                
                                                
                                                
        