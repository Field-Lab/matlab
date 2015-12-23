

%blueberry
path_and_name{1,1} = '/Volumes/lab/Experiments/Array/Analysis/Chichilnisky-lab/2011-09-23-0/data004/data004';
path_and_name{1,2} = '2011-09-23-0_rf-4-bayes-msf_10.00-10-0001-rat';


% load data
datarun = load_data(path_and_name{1,1});
datarun = load_params(datarun,struct('verbose',1));  
datarun = load_sta(datarun,'load_sta',[]);
datarun = set_polarities(datarun);
datarun = import_single_cone_data(datarun, path_and_name{1,2});    
datarun.cones.mosaic = make_mosaic_struct(datarun.cones.centers);

datarun = get_sta_summaries(datarun, 'all')

datarun.stimulus.x_start = 0;
datarun.stimulus.x_end = 200;
datarun.stimulus.y_start = 0;
datrun.stimulus.y_end = 200;
datarun.field_height = 200;
datarun.field_width = 200;


cell_spec = [346 1818 3934 5656 6906 7415 7756]; % type 2
cell_spec = [466 4370 4457 5479 7456]; % type 4
cell_spec = [949 2781 4221 4279 5273 6273 6275 6621]; % type 1

temp_polarity = 1;

% plot background

% get size and color
%y = datarun.stimulus.field_height;
%x = datarun.stimulus.field_width;
y = 200;
x = 200;
rgb = [.35 .35 .35];

% generate background matrix
plot_mat = cat(3,repmat(rgb(1),y,x),repmat(rgb(2),y,x),repmat(rgb(3),y,x));

% plot it
figure(10);clf;image(plot_mat);axis image; hold on


clear selection_params
selection_params.thresh = 0.2;
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
    'thresh', selection_params.thresh)


% plot center cone halos
cone_colors = [0.5 0.5 0.5] + (0.5*temp_polarity); 
cone_colors = repmat(cone_colors,4,1);

plot_cone_mosaic(datarun,'fig_or_axes',10,'cone_size',8,'clear',0,'bg_color',[],...
                'cone_roi', cone_center_roi, 'roi_highlight', false, 'type_colors', cone_colors)

% plot center cones
plot_cone_mosaic(datarun,'fig_or_axes',10,'cone_size',6,'clear',0,'bg_color',[])



print(10,'~/Desktop/spider','-dpdf')


%%  plot RF fits for the population
datarun = get_sta_fits_from_vision(datarun, 'all')

plot_rf_summaries(datarun, cell_spec, 'foa', 1, 'plot_fits', true)
print(1,'~/Desktop/rf-outlines-all','-dpdf')

%% plot single rf, best example from each type

%small
plot_rf(datarun, 6275, 'scale', 8)
print(1,'~/Desktop/example-small','-dpdf')

%medium
plot_rf(datarun, 7456, 'scale', 8)
print(1,'~/Desktop/example-medium','-dpdf')

%large
plot_rf(datarun, 7415, 'scale', 8)
print(1,'~/Desktop/example-large','-dpdf')


%%

temp_rf = get_rf(datarun, 7428);
rf_filter = make_Gaussian_two_d('sd_x', 0.75, 'sd_y', 0.75);

filt_rf = conv2(temp_rf, rf_filter, 'same');

image(norm_image(filt_rf))





