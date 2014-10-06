% 
% master_data_path = '/snle/lab/Experiments/Array/Analysis/2008-08-27-5/data003/data003/data003';
% slave_data_path = '/snle/lab/Experiments/Array/Analysis/2008-08-27-5/data005/data005/data005';

master_data_path = '/Volumes/Analysis/2008-08-27-5/data003/data003/data003';
slave_data_path = '/Volumes/Analysis/2008-08-27-5/data005/data005/data005';


array_type = 519;

% load master data
clear temp_datarun
temp_datarun = load_data(master_data_path);
temp_datarun = load_params(temp_datarun,struct('verbose',1));  
temp_datarun = load_ei(temp_datarun, 'all', struct('array_type', array_type));
datarun{1} = temp_datarun;

% load slave data
clear temp_datarun
temp_datarun = load_data(slave_data_path);
temp_datarun = load_neurons(temp_datarun);
temp_datarun = load_ei(temp_datarun, 'all',struct('array_type', array_type));
datarun{2} = temp_datarun;


threshold = 0.92; % correlaton threhold for mapping
cell_type = {4};

% EI mapping from RF run to control condition
master_cell_type = cell_type;
slave_cell_type = 'all';
[datarunA, datarunB] = map_gdf(datarun{1}, datarun{2}, 'corr_threshold', threshold, 'master_cell_type', master_cell_type, 'slave_cell_type', slave_cell_type, 'verbose', true, 'troubleshoot', true);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% plantain
% path_and_name{1,1} = '/snle/lab/Experiments/Array/Analysis/2008-08-27-5/data003/data003/data003';
% path_and_name{1,2} = '2008-08-27-5_data003_data003_data003-bayes-msf_70.00';

path_and_name{1,1} = '/Volumes/Analysis/2008-08-27-5/data003/data003/data003';
path_and_name{1,2} = '2008-08-27-5_data003_data003_data003-bayes-msf_70.00--standard';

% load data
datarun = load_data(path_and_name{1,1});
datarun = load_params(datarun,struct('verbose',1));  
datarun = load_sta(datarun,'load_sta',[]);
datarun = set_polarities(datarun);
datarun = import_single_cone_data(datarun, path_and_name{1,2});    
datarun.cones.mosaic = make_mosaic_struct(datarun.cones.centers);


cell_id = 1381;
cell_index = get_cell_indices(datarun, cell_id);
datarun = get_sta_summaries(datarun, cell_id);


figure(1); clf;
plot_rf(datarun, cell_id, 'scale', 20, 'polarity', true, 'foa', 1)
plot_cone_mosaic(datarun , 'cone_size', 3, 'fig_or_axes', 1, 'clear', false, 'bg_color', [])
print(1, '~/Desktop/rf.pdf', '-dpdf')


figure(2); clf;
temp_rf = datarun.stas.rfs{cell_index};
temp_rf(temp_rf > 0) = 0.0;
imagesc(norm_image(-1*temp_rf))

rect_datarun = datarun;
rect_datarun.stas.rfs{cell_index} = temp_rf;

plot_rf(rect_datarun, cell_id, 'scale', 10, 'polarity', true, 'foa', 2)
plot_cone_mosaic(datarun , 'cone_size', 3, 'fig_or_axes', 2, 'clear', false, 'bg_color', [])
print(2, '~/Desktop/surround.pdf', '-dpdf')





%%%% select cell

% OFF cells
%cell_id = 1099; % M off, L ON
%cell_id = 3136; % L OFF, M, OFF
% On cells
cell_id = 813; % L ON, M OFF

%cell_id = 1381 % cell of figure 2h


% other cells
%cell_id =  976;
%cell_id = 1666;
%cell_id = 571;
%cell_id = 766;
%cell_id = 1652;  % M ON
%cell_id = 1411;  % nice L off cell


% non opponent
%cell_id = 3440


%cell_id = datarun.cell_types{4}.cell_ids(83)

cell_index = get_cell_indices(datarun,cell_id);

cell_threshold = 0.033;

% plot background
% get size and color
y = datarun.stimulus.field_height;
x = datarun.stimulus.field_width;
rgb = [0.35 0.35 0.35];
% generate background matrix
plot_mat = cat(3,repmat(rgb(1),y,x),repmat(rgb(2),y,x),repmat(rgb(3),y,x));
% plot it
figure(10);clf;image(plot_mat);axis image; hold on


% plot_center and surround spiders
plot_cell_sampling(datarun,cell_id,'type','spider','fig_or_axes',10,'clear',0,'plot_cones',0,...
   'line_width',[realmin 2.5],'cell_colors',[1 1 1],'plot_radius',[0 6],'thresh',cell_threshold,'polarity',0,'contiguity',true)


type_colors = [1 0 0; 0 1 0; 0 0 1; 0 0 0];

% plot all cones
plot_cone_mosaic(datarun,'fig_or_axes',10,'cone_size',4,'clear',0,'bg_color',[], 'type_colors', type_colors)





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get center information and plot center cones
% get center cones
clear selection_params
selection_params.radius = [0 6];
selection_params.contiguity = true;
selection_params.thresh = 0.03;
selection_params.polarity = 1;
[weights, center_selection, extras] = select_cone_weights(datarun, cell_id, selection_params);

% plot center cone halos
%cone_colors = [0.5 0.5 0.5] + (0.5 * repmat(datarun.stas.polarities{cell_index}, 1,3));
cone_colors = [1 1 1];
cone_colors = repmat(cone_colors,4,1);

plot_cone_mosaic(datarun,'fig_or_axes',10,'cone_size',7,'clear',0,'bg_color',[],...
                'cone_roi', center_selection, 'roi_highlight', false)

% plot center cones
plot_cone_mosaic(datarun,'fig_or_axes',10,'cone_size',5,'clear',0,'bg_color',[],...
                'cone_roi', center_selection, 'roi_highlight', false)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get center information and plot center cones
% get center cones
clear selection_params
selection_params.radius = [0 6];
selection_params.contiguity = true;
selection_params.thresh = cell_threshold;
selection_params.polarity = -1;
[weights, center_selection, extras] = select_cone_weights(datarun, cell_id, selection_params);

% plot surround cone halos
%cone_colors = [0.5 0.5 0.5] + (0.5 * repmat(datarun.stas.polarities{cell_index}, 1,3));
cone_colors = [0 0 0];
cone_colors = repmat(cone_colors,4,1);

plot_cone_mosaic(datarun,'fig_or_axes',10,'cone_size',7,'clear',0,'bg_color',[],...
                'cone_roi', center_selection, 'roi_highlight', false, 'type_colors', cone_colors)

% plot center cones
plot_cone_mosaic(datarun,'fig_or_axes',10,'cone_size',5,'clear',0,'bg_color',[], 'cone_roi', center_selection, 'roi_highlight', false)

print(10,'/snle/home/gfield/Desktop/spider','-dpdf')







