% plantain
path_and_name{1,1} = '/snle/lab/Experiments/Array/Analysis/2008-08-27-5/data003/data003/data003';
path_and_name{1,2} = '2008-08-27-5_data003_data003_data003-bayes-msf_70.00--standard';
%path_and_name{1,2} = 'erroneous_normalization/2008-08-27-5_data003_data003_data003-bayes-msf_70.00--standard-old';


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
    'thresh', selection_params.thresh,...
     'label', true)


% plot center cone halos
cone_colors = [0.5 0.5 0.5] + (0.5*temp_polarity); 
cone_colors = repmat(cone_colors,4,1);

plot_cone_mosaic(datarun,'fig_or_axes',10,'cone_size',12,'clear',0,'bg_color',[],...
                'cone_roi', cone_center_roi, 'roi_highlight', false, 'type_colors', cone_colors)

% plot center cones
plot_cone_mosaic(datarun,'fig_or_axes',10,'cone_size',6,'clear',0,'bg_color',[])


%%
% cells of interest
window_size = 15;

off_mid_one = 1412;
off_mid_two = 1411;
on_midget_one = 1456;

% example cell 1
cell_one = on_midget_one;
cell_one_index = get_cell_indices(datarun, cell_one);
plot_rf(datarun, cell_one, 'foa', 1, 'scale', 8)

%get COM
temp_COM = datarun.stas.rf_coms{cell_one_index};
x_min = temp_COM(1) - window_size;
x_max = temp_COM(1) + window_size;
y_min = temp_COM(2) - window_size;
y_max = temp_COM(2) + window_size;
axis([x_min, x_max, y_min, y_max])

%%
on_midget_two = 931;
% example cell 1
cell_one = on_midget_two;
cell_one_index = get_cell_indices(datarun, cell_one);
plot_rf(datarun, cell_one, 'foa', 2, 'scale', 8)

%get COM
temp_COM = datarun.stas.rf_coms{cell_one_index};
x_min = temp_COM(1) - window_size;
x_max = temp_COM(1) + window_size;
y_min = temp_COM(2) - window_size;
y_max = temp_COM(2) + window_size;
axis([x_min, x_max, y_min, y_max])

%%
on_midget_two = 1276;
% example cell 1
cell_one = on_midget_two;
cell_one_index = get_cell_indices(datarun, cell_one);
plot_rf(datarun, cell_one, 'foa', 3, 'scale', 8)

%get COM
temp_COM = datarun.stas.rf_coms{cell_one_index};
x_min = temp_COM(1) - window_size;
x_max = temp_COM(1) + window_size;
y_min = temp_COM(2) - window_size;
y_max = temp_COM(2) + window_size;
axis([x_min, x_max, y_min, y_max])
%%
off_midget_one = 1742;
% example cell 1
cell_one = off_midget_one;
cell_one_index = get_cell_indices(datarun, cell_one);
plot_rf(datarun, cell_one, 'foa', 7, 'scale', 8)

%get COM
temp_COM = datarun.stas.rf_coms{cell_one_index};
x_min = temp_COM(1) - window_size;
x_max = temp_COM(1) + window_size;
y_min = temp_COM(2) - window_size;
y_max = temp_COM(2) + window_size;
axis([x_min, x_max, y_min, y_max])

%%
off_midget_two = 1411;
% example cell 1
cell_one = off_midget_two;
cell_one_index = get_cell_indices(datarun, cell_one);
plot_rf(datarun, cell_one, 'foa', 5, 'scale', 8)

%get COM
temp_COM = datarun.stas.rf_coms{cell_one_index};
x_min = temp_COM(1) - window_size;
x_max = temp_COM(1) + window_size;
y_min = temp_COM(2) - window_size;
y_max = temp_COM(2) + window_size;
axis([x_min, x_max, y_min, y_max])
%%
off_parasol_one = 1686;
window_size = 20;
% example cell 1
cell_one = off_parasol_one;
cell_one_index = get_cell_indices(datarun, cell_one);
plot_rf(datarun, cell_one, 'foa', 6, 'scale', 8)
plot_cone_mosaic(datarun, 'fig_or_axes', 6, 'cone_size', 6, 'bg_color', [], 'clear', false)


%get COM
temp_COM = datarun.stas.rf_coms{cell_one_index};
x_min = temp_COM(1) - window_size;
x_max = temp_COM(1) + window_size;
y_min = temp_COM(2) - window_size;
y_max = temp_COM(2) + window_size;
axis([x_min, x_max, y_min, y_max])

%%
% kiwi
path_and_name{1,1} = '/snle/lab/Experiments/Array/Analysis/2008-05-13-3/data006/data006';
path_and_name{1,2} = '2008-05-13-3_data006_data006-bayes-msf_85.00--standard';
%path_and_name{1,2} = 'erroneous_normalization/2008-08-27-5_data003_data003_data003-bayes-msf_70.00--standard-old';


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

fig_num = 11;

% generate background matrix
plot_mat = cat(3,repmat(rgb(1),y,x),repmat(rgb(2),y,x),repmat(rgb(3),y,x));

% plot it
figure(fig_num);clf;image(plot_mat);axis image; hold on


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
plot_cell_sampling(datarun,cell_spec,'type','spider','fig_or_axes',fig_num,'clear',0,'plot_cones',0,...
   'line_width',[realmin 2.0],'cell_colors',...
    [0.5 0.5 0.5] + (0.5*temp_polarity),'plot_radius',[0 6],...
    'thresh', selection_params.thresh,...
     'label', true)


% plot center cone halos
cone_colors = [0.5 0.5 0.5] + (0.5*temp_polarity); 
cone_colors = repmat(cone_colors,4,1);

plot_cone_mosaic(datarun,'fig_or_axes',fig_num,'cone_size',12,'clear',0,'bg_color',[],...
                'cone_roi', cone_center_roi, 'roi_highlight', false, 'type_colors', cone_colors)

% plot center cones
plot_cone_mosaic(datarun,'fig_or_axes',fig_num,'cone_size',6,'clear',0,'bg_color',[])

