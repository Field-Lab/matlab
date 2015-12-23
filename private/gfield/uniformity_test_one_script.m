LMS_mean = [0 0 0];
LMS_cov = [1 0.8 0.4; 0.8 1 0.5; 0.4 0.5 1];
num_LMS_samples = 100;

rand_LMS_triplets = mvnrnd(LMS_mean, LMS_cov, num_LMS_samples);

%% load a dataset and cone info

% plantain
path_and_name{1,1} = '/snle/lab/Experiments/Array/Analysis/2008-08-27-5/data003/data003/data003';
path_and_name{1,2} = '2008-08-27-5_data003_data003_data003-bayes-msf_70.00--standard';

clear datarun 
datarun = load_data(path_and_name{1,1});
datarun = load_params(datarun,struct('verbose',1));  
datarun = load_sta(datarun,'load_sta',[]);
datarun = import_single_cone_data(datarun, path_and_name{1,2});    
datarun.cones.mosaic = make_mosaic_struct(datarun.cones.centers);

cell_type = {3};

temp_cell_indices = get_cell_indices(datarun, cell_type);

[mosaic_weights, selection, extras] = select_cone_weights(datarun, cell_type,...
                                            'thresh', 0.05,...
                                            'radius', [0 inf], 'polarity', 1,...
                                            'contiguity', true,'scale', 3.0,...
                                            'remove_cones', 'S');   

num_rgcs = length(temp_cell_indices);
new_datarun = extras.new_datarun;

L_indices = find(new_datarun.cones.types == 'L');
M_indices = find(new_datarun.cones.types == 'M');

total_signal = zeros(num_LMS_samples, num_rgcs);
for rgc = 1:num_rgcs
    temp_weights = mosaic_weights(:,rgc) .* selection(:,rgc);
    temp_weights = temp_weights ./ sum(temp_weights);
    temp_l_total = sum(temp_weights(L_indices)) * rand_LMS_triplets(:,1);
    temp_m_total = sum(temp_weights(M_indices)) * rand_LMS_triplets(:,2);
    
    total_signal(:,rgc) = temp_l_total + temp_m_total;
end
                                            
var_signals = var(total_signal,0, 1);    

mean_var_signals = mean(var_signals)

%% permute cones and recompute
num_iters = 100;
cone_type_list = new_datarun.cones.types;
num_cones = length(cone_type_list);
permuted_cone_var_signals = zeros(num_iters, num_rgcs);

for iter = 1:num_iters
    temp_type_list = cone_type_list(randperm(num_cones));
    new_L_indices = find(temp_type_list == 'L');
    new_M_indices = find(temp_type_list == 'M');
    
    total_signal = zeros(num_LMS_samples, num_rgcs);
    for rgc = 1:num_rgcs
        temp_weights = mosaic_weights(:,rgc) .* selection(:,rgc);
        temp_weights = temp_weights ./ sum(temp_weights);
        temp_l_total = sum(temp_weights(new_L_indices)) * rand_LMS_triplets(:,1);
        temp_m_total = sum(temp_weights(new_M_indices)) * rand_LMS_triplets(:,2);

        total_signal(:,rgc) = temp_l_total + temp_m_total;
    end    

    permuted_cone_var_signals(iter,:) = var(total_signal,0,1);  
    
end

mean_var_permuted_signals = mean(permuted_cone_var_signals, 2);

mean_over_iters = mean(mean_var_permuted_signals)



    
    