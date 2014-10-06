% analyze the convergence of cones across multiple RGC types to identify
% recordings with large (small) cone convergence

clear 

% set analysis parameters
binarize_weights_flag = false; % if true, weights are binarized
use_cone_files = true; % if true, cone files are loaded from server
clumped_flag = true; % if true, clumped cone mosaics are loaded, 'use_cone_files', must also be true
rf_center_threshold = 0.1; % defines threshold for what is RF center
cell_cutoff = 50; % number of rgcs that are required to analyze the mosaic
num_permuted_cone_mosaics = 100;  % num
verbose = 0;


% get paths to data of interest and number of datasets
[LMS_paths, LMS_names, cone_paths] = get_LMS_paths('high', 'cone_finding', 'standard');
num_datasets = length(LMS_paths);

cell_types = [1,2,3,4];
num_cell_types = length(cell_types);

cone_convergences = zeros(num_datasets, num_cell_types);
for dataset = 1:num_datasets

    % load information from a data set
    clear datarun
    datarun = load_data(LMS_paths{dataset});
    datarun = load_params(datarun,struct('verbose',1));  
    datarun = load_sta(datarun,'load_sta',[]);
    datarun = import_single_cone_data(datarun, cone_paths{dataset}, 'overwrite_cell_types', true);    
    datarun.cones.mosaic = make_mosaic_struct(datarun.cones.centers);
    
    for ctype = 1:length(cell_types)
        
        cell_type = cell_types(ctype);

        % extract connectivity
        [mosaic_weights, selection, extras] = select_cone_weights(datarun, {cell_type},...
                                                    'thresh', rf_center_threshold,...
                                                    'radius', [0 inf], 'polarity', 1,...
                                                    'contiguity', true,'scale', 3.0,...
                                                    'remove_cones', 'S');   
        connectivity = mosaic_weights .* selection; % keep weights continuous valued

        [num_cones, num_rgcs] = size(connectivity);
        
        temp_cone_conv = zeros(1,num_rgcs);
        for rgc = 1:num_rgcs
            temp_cone_conv(rgc) = length(find(connectivity(:,rgc)));
        end
        
        cone_convergences(dataset, ctype) = median(temp_cone_conv);
        
    end
end
        
        
cone_convergences


% plantain
path_and_name{1,1} = '/snle/lab/Experiments/Array/Analysis/2008-08-27-5/data003/data003/data003';
path_and_name{1,2} = '2008-08-27-5_data003_data003_data003-bayes-msf_70.00--standard';
%path_and_name{1,2} = 'erroneous_normalization/2008-08-27-5_data003_data003_data003-bayes-msf_70.00--standard-old';

%  plantain alternative S cone threshold
path_and_name{1,1} = '/snle/lab/Experiments/Array/Analysis/2008-08-27-5/data003/data003/data003';
path_and_name{1,2} = '2008-08-27-5_data003_data003_data003-bayes-msf_LM_70.00-msf_S_20.00---more_blues';

%  plantain old cone finding
path_and_name{1,1} = '/snle/lab/Experiments/Array/Analysis/2008-08-27-5/data003/data003/data003';
path_and_name{1,2} = 'plantain-old';

%  plantain old cone finding
path_and_name{1,1} = '/snle/lab/Experiments/Array/Analysis/2008-08-27-5/data003/data003/data003';
path_and_name{1,2} = 'plantain-old';

%blueberry
path_and_name{1,1} = '/snle/lab/Experiments/Array/Analysis/2008-08-26-2/data001-s6369-s9551/data001-s6369-s9551';
path_and_name{1,2} = '2008-08-26-2_data001-s6369-s9551_data001-s6369-s9551-bayes-msf_20.00--standard';


% grapes
path_and_name{1,1} = '/snle/lab/Experiments/Array/Analysis/2007-03-27-2/data014/data014/data014';
path_and_name{1,2} = '2007-03-27-2_data014_data014_data014-bayes-msf_25.00--standard';


% load data
datarun = load_data(path_and_name{1,1});
datarun = load_params(datarun,struct('verbose',1));  
datarun = load_sta(datarun,'load_sta',[]);
datarun = set_polarities(datarun);
datarun = import_single_cone_data(datarun, path_and_name{1,2});    
datarun.cones.mosaic = make_mosaic_struct(datarun.cones.centers);


cell_spec = {4};
temp_polarity = 1;

% plot background

% get size and color
y = datarun.stimulus.field_height;
x = datarun.stimulus.field_width;
rgb = [.35 .35 .35];

% generate background matrix
plot_mat = cat(3,repmat(rgb(1),y,x),repmat(rgb(2),y,x),repmat(rgb(3),y,x));

% plot it
figure(10);clf;image(plot_mat);axis image; hold on


clear selection_params
selection_params.thresh = 0.1;
selection_params.radius = [0 inf];
selection_params.contiguity = true;
selection_params.polarity = 1;
[weights, center_selection, extras] = select_cone_weights(datarun, cell_spec, selection_params);

summed_selection = sum(center_selection, 2);
sampled_indices = find(summed_selection > 0);
cone_center_roi = zeros(length(summed_selection),1);
cone_center_roi(sampled_indices) = 1;

% plot spiders
plot_cell_sampling(datarun,cell_spec,'type','spider','fig_or_axes',10,'clear',0,'plot_cones',0,...
   'line_width',[realmin 2.0],'cell_colors',...
    [0.5 0.5 0.5] + (0.5*temp_polarity),'plot_radius',[0 6],...
    'thresh', selection_params.thresh, 'label', true)


% plot center cone halos
cone_colors = [0.5 0.5 0.5] + (0.5*temp_polarity); 
cone_colors = repmat(cone_colors,4,1);

plot_cone_mosaic(datarun,'fig_or_axes',10,'cone_size',3,'clear',0,'bg_color',[],...
                'cone_roi', cone_center_roi, 'roi_highlight', false, 'type_colors', cone_colors)

% plot center cones
plot_cone_mosaic(datarun,'fig_or_axes',10,'cone_size',4,'clear',0,'bg_color',[])



print(10,'~/Desktop/spider','-dpdf')
