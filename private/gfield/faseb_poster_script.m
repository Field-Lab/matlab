path_and_name{1,1} = '/snle/lab/Experiments/Array/Analysis/2008-08-27-5/data003/data003/data003';
path_and_name{1,2} = '2008-08-27-5_data003_data003_data003-bayes-msf_70.00--standard';

path_and_name{1,1} = '/snle/lab/Experiments/Array/Analysis/2010-03-05-2/data013/data013';
path_and_name{1,2} = '2010-03-05-2_rf-13-apple-gf-bayes-msf_25.00--standard';

% load data
datarun = load_data(path_and_name{1,1});
datarun = load_params(datarun,struct('verbose',1));  
datarun = load_sta(datarun,'load_sta',[]);
datarun = set_polarities(datarun);
datarun = import_single_cone_data(datarun, path_and_name{1,2});    
datarun.cones.mosaic = make_mosaic_struct(datarun.cones.centers);


datarun = get_sta_summaries(datarun, {1,2,3,4}, 'keep_stas', false);


plot_rf_portraits(datarun, 947, 'plot_radius', 25, 'figure', 1, 'scale_factor', 10)
print(1, '~/Desktop/raw-sta.pdf', '-dpdf')


cell_spec = 947;
temp_polarity = 1;

% plot background

% get size and color
y = datarun.stimulus.field_height;
x = datarun.stimulus.field_width;
rgb = [0.5 0.5 0.5];

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
   'line_width',[realmin 4.0],'cell_colors',...
    [0.5 0.5 0.5] + (0.5*temp_polarity),'plot_radius',[0 6],...
    'thresh', selection_params.thresh)


% plot center cone halos
cone_colors = [0.5 0.5 0.5] + (0.5*temp_polarity); 
cone_colors = repmat(cone_colors,4,1);

plot_cone_mosaic(datarun,'fig_or_axes',10,'cone_size',5,'clear',0,'bg_color',[],...
                'cone_roi', cone_center_roi, 'roi_highlight', false, 'type_colors', cone_colors)

% plot center cones
plot_cone_mosaic(datarun,'fig_or_axes',10,'cone_size',4,'clear',0,'bg_color',[])

print(10,'~/Desktop/spider','-dpdf')



