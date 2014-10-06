% plantain
path_and_name{1,1} = '/snle/lab/Experiments/Array/Analysis/2008-08-27-5/data003/data003/data003';
path_and_name{1,2} = '2008-08-27-5_data003_data003_data003-bayes-msf_70.00--standard';


% load data
datarun = load_data(path_and_name{1,1});
datarun = load_params(datarun,struct('verbose',1));  
datarun = load_sta(datarun,'load_sta',[]);
datarun = set_polarities(datarun);
datarun = import_single_cone_data(datarun, path_and_name{1,2});    
datarun.cones.mosaic = make_mosaic_struct(datarun.cones.centers);



cell_id = 1456; % on midget cell
%cell_id = 1336; % off midget cell
%cell_id = 1412; % off midget cell

cell_index = get_cell_indices(datarun,cell_id);

cell_threshold = 0.1;

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
   'line_width',[realmin 2.5],'cell_colors',[1 1 1],'plot_radius',[0 6],'thresh',cell_threshold,'polarity',1,'contiguity',true)

de_sta = 0;
type_colors = [1 de_sta de_sta; de_sta 1 de_sta; 0 0 1; 0 0 0];

% plot all cones
plot_cone_mosaic(datarun,'fig_or_axes',10,'cone_size',4,'clear',0,'bg_color',[], 'type_colors', type_colors)


% get center information and plot center cones
% get center cones
clear selection_params
selection_params.radius = [0 6];
selection_params.contiguity = false;
selection_params.thresh = cell_threshold;
selection_params.polarity = 1;
[weights, center_selection, extras] = select_cone_weights(datarun, cell_id, selection_params);

% plot center cone halos
%cone_colors = [0.5 0.5 0.5] + (0.5 * repmat(datarun.stas.polarities{cell_index}, 1,3));
cone_colors = [1 1 1];
cone_colors = repmat(cone_colors,4,1);

plot_cone_mosaic(datarun,'fig_or_axes',10,'cone_size',7,'clear',0,'bg_color',[],...
                'cone_roi', center_selection, 'roi_highlight', false, 'type_colors', cone_colors)

% plot center cones
plot_cone_mosaic(datarun,'fig_or_axes',10,'cone_size',5,'clear',0,'bg_color',[],...
                'cone_roi', center_selection, 'roi_highlight', false, 'type_colors', type_colors)


print(10,'/snle/home/gfield/Desktop/spider','-dpdf')




