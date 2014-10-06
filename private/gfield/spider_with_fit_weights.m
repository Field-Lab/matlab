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



datarun = fit_cone_rfs(datarun,{1,2,3,4});


cell_types = {4};
rgc_indices = get_cell_indices(datarun, cell_types);
num_rgcs = length(rgc_indices);
num_cones = length(datarun.cones.types);


% get weights from RF fits
new_mosaic_weights = zeros(num_cones, num_rgcs);

old_datarun = datarun;

for rgc = 1:num_rgcs
    fit_params = [];
    fit_params.center = datarun.cones.rf_fits{rgc_indices(rgc)}.center;
    fit_params.center_scale = datarun.cones.rf_fits{rgc_indices(rgc)}.center_scale;
    fit_params.center_radius = datarun.cones.rf_fits{rgc_indices(rgc)}.center_radius;
    fit_params.surround_scale = datarun.cones.rf_fits{rgc_indices(rgc)}.surround_scale;
    fit_params.surround_radius = datarun.cones.rf_fits{rgc_indices(rgc)}.surround_radius;
    
    temp_weights = dog_fit(datarun.cones.centers, fit_params);

    new_mosaic_weights(:, rgc) = temp_weights;
    
    datarun.cones.weights(:, rgc_indices(rgc)) = temp_weights;
end

datarun = new_datarun;
figure_num = 10;

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
figure(figure_num);clf;image(plot_mat);axis image; hold on


clear selection_params
selection_params.thresh = 0.05;
selection_params.radius = [0 inf];
selection_params.contiguity = true;
selection_params.polarity = 1;
[weights, center_selection, extras] = select_cone_weights(datarun, cell_spec, selection_params);

summed_selection = sum(center_selection, 2);
sampled_indices = find(summed_selection > 0);
cone_center_roi = zeros(length(summed_selection),1);
cone_center_roi(sampled_indices) = 1;

% plot spiders
plot_cell_sampling(datarun,cell_spec,'type','spider','fig_or_axes',figure_num,'clear',0,'plot_cones',0,...
   'line_width',[realmin 2.0],'cell_colors',...
    [0.5 0.5 0.5] + (0.5*temp_polarity),'plot_radius',[0 6],...
    'thresh', selection_params.thresh)


% plot center cone halos
cone_colors = [0.5 0.5 0.5] + (0.5*temp_polarity); 
cone_colors = repmat(cone_colors,4,1);

plot_cone_mosaic(datarun,'fig_or_axes',figure_num,'cone_size',4,'clear',0,'bg_color',[],...
                'cone_roi', cone_center_roi, 'roi_highlight', false, 'type_colors', cone_colors)

% plot center cones
plot_cone_mosaic(datarun,'fig_or_axes',figure_num,'cone_size',3,'clear',0,'bg_color',[])

print(10,'~/Desktop/spider','-dpdf')

%%
%load fits from coarse stixel run
data_path = '/snle/lab/Experiments/Array/Analysis/2008-08-27-5/data008/data008/data008';
obvius_fit_path = '/snle/lab/Experiments/Array/Analysis/2008-08-27-5/rf-8-gf/';

datarun_fit = load_data(data_path);
datarun_fit = load_params(datarun_fit);
datarun_fit = load_neurons(datarun_fit);
datarun_fit = load_index(datarun_fit);
datarun_fit.names.obvius_fit_path = obvius_fit_path;

datarun_fit = load_obvius_sta_fits(datarun_fit);
datarun_fit = get_sta_fits_from_obvius(datarun_fit, {1,2,3,4,5});


cell_types = {1,2,3,4};
pixel_rescale = 5; % scale factor between 1x1 and coarse spatial STA run

temp_indices = get_cell_indices(datarun_fit, cell_types);
num_rgcs = length(temp_indices);

cone_locations = datarun.cones.centers;
new_datarun = datarun;
new_weights = zeros(num_cones,num_rgcs);
for rgc = 1:num_rgcs
    % extract and organized fit information
    the_fit = datarun_fit.obvius.sta_fits{temp_indices(rgc)};
    fit_params = [];
    fit_params.center = (the_fit.mean + [1,1]) * pixel_rescale;
    fit_params.center_scale = the_fit.scale(1);
    fit_params.center_sd = the_fit.sd .* pixel_rescale;
    fit_params.angle = the_fit.angle;
    fit_params.surround_sd_scale = the_fit.surround_sd_scale;
    fit_params.surround_scale = the_fit.surround_scale;
    
    temp_weights = two_d_dog_fit(cone_locations, fit_params);
    new_weights(:,rgc) = temp_weights;
    new_datarun.cones.weights(:,temp_indices(rgc)) = temp_weights;
end

    
figure_num = 10;

cell_spec = {3};
temp_polarity = 1;

% plot background

% get size and color
y = datarun.stimulus.field_height;
x = datarun.stimulus.field_width;
rgb = [.35 .35 .35];

% generate background matrix
plot_mat = cat(3,repmat(rgb(1),y,x),repmat(rgb(2),y,x),repmat(rgb(3),y,x));

% plot it
figure(figure_num);clf;image(plot_mat);axis image; hold on


clear selection_params
selection_params.thresh = 0.01;
selection_params.radius = [0 inf];
selection_params.contiguity = false;
selection_params.polarity = 1;
[weights, center_selection, extras] = select_cone_weights(datarun, cell_spec, selection_params);

summed_selection = sum(center_selection, 2);
sampled_indices = find(summed_selection > 0);
cone_center_roi = zeros(length(summed_selection),1);
cone_center_roi(sampled_indices) = 1;

% plot spiders
plot_cell_sampling(new_datarun,cell_spec,'type','spider','fig_or_axes',figure_num,'clear',0,'plot_cones',0,...
   'line_width',[realmin 2.0],'cell_colors',...
    [0.5 0.5 0.5] + (0.5*temp_polarity),'plot_radius',[0 6],...
    'thresh', selection_params.thresh,'contiguity', selection_params.contiguity)


% plot center cone halos
cone_colors = [0.5 0.5 0.5] + (0.5*temp_polarity); 
cone_colors = repmat(cone_colors,4,1);

plot_cone_mosaic(datarun,'fig_or_axes',figure_num,'cone_size',4,'clear',0,'bg_color',[],...
                'cone_roi', cone_center_roi, 'roi_highlight', false, 'type_colors', cone_colors)

% plot center cones
plot_cone_mosaic(datarun,'fig_or_axes',figure_num,'cone_size',8,'clear',0,'bg_color',[])

    





















