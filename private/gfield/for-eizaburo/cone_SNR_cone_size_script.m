% get paths to data of interest and number of datasets
% blueberry
data_path = '/snle/lab/Experiments/Array/Analysis/2008-08-26-2/data001-s6369-s9551/data001-s6369-s9551';
cone_path{1} = '2008-08-26-2_rf-1-blueberry-bayes-msf_20.00-.70xcone_radius';
cone_path{2} = '2008-08-26-2_rf-1-blueberry-bayes-msf_20.00-.85xcone_radius';
cone_path{3} = '2008-08-26-2_data001-s6369-s9551_data001-s6369-s9551-bayes-msf_20.00--standard';
cone_path{4} = '2008-08-26-2_rf-1-blueberry-bayes-msf_20.00-1.5xcone_radius';
cone_path{5} = '2008-08-26-2_rf-1-blueberry-bayes-msf_20.00-2.0xcone_radius';

% apple
data_path = '/snle/lab/Experiments/Array/Analysis/2010-03-05-2/data013/data013';
cone_path{1} = '2010-03-05-2_rf-13-apple-gf-bayes-msf_25.00-.70xcone_radius';
cone_path{2} = '2010-03-05-2_rf-13-apple-gf-bayes-msf_25.00-.85xcone_radius';
cone_path{3} = '2010-03-05-2_rf-13-apple-gf-bayes-msf_25.00--standard';
cone_path{4} = '2010-03-05-2_rf-13-apple-gf-bayes-msf_25.00-1.5xcone_radius';
cone_path{5} = '2010-03-05-2_rf-13-apple-gf-bayes-msf_25.00-2.0xcone_radius';

% plantain
data_path = '/snle/lab/Experiments/Array/Analysis/2008-08-27-5/data003/data003/data003';
cone_path{1} = '2008-08-27-5_data003_data003_data003-bayes-msf_70.00-0.70xcone_radius';
cone_path{2} = '2008-08-27-5_data003_data003_data003-bayes-msf_70.00-0.85xcone_radius';
cone_path{3} = '2008-08-27-5_data003_data003_data003-bayes-msf_70.00--standard';
cone_path{4} = '2008-08-27-5_data003_data003_data003-bayes-msf_70.00-1.5xcone_radius';
cone_path{5} = '2008-08-27-5_data003_data003_data003-bayes-msf_70.00-2.0xcone_radius';

num_datasets = length(cone_path);

mosaic_counter = 0;
rf_center_threshold = 0.1;
dataset_SNRs = zeros(1,num_datasets);

for dataset = 1:num_datasets

    % load information from a data set
    clear datarun
    datarun = load_data(data_path);
    datarun = load_params(datarun,struct('verbose',1));  
    datarun = load_sta(datarun,'load_sta',[]);
    datarun = import_single_cone_data(datarun, cone_path{dataset}, 'overwrite_cell_types', true);    
    datarun.cones.mosaic = make_mosaic_struct(datarun.cones.centers);
    


    % extract connectivity
    [mosaic_weights, selection, extras] = select_cone_weights(datarun, {3,4},...
                                                'thresh', rf_center_threshold,...
                                                'radius', [0 inf], 'polarity', 1,...
                                                'contiguity', true,'scale', 3.0,...
                                                'remove_cones', 'S');   


    temp_SNRs = get_cone_weight_SNR(mosaic_weights, selection);  
    dataset_SNRs(dataset) = mean(temp_SNRs);
    
end

dataset_SNRs