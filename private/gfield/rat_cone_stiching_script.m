%% load data
% rat data
path_and_name{1,1} = '/snle/lab/Experiments/Array/Analysis/2011-09-23-0/data004/data004';
path_and_name{1,2} = '2011-09-23-0_rf-4-bayes-msf_10.00-10-0001-rat';


% load data
datarun = load_data(path_and_name{1,1});
datarun = load_params(datarun,struct('verbose',1));  
datarun = load_sta(datarun,'load_sta',[]);
datarun = set_polarities(datarun);
datarun = import_single_cone_data(datarun, path_and_name{1,2});    
datarun.cones.mosaic = make_mosaic_struct(datarun.cones.centers);

datarun = get_sta_summaries(datarun, {8,9});

datarun.stimulus.x_start = 0;
datarun.stimulus.x_end = 200;
datarun.stimulus.y_start = 0;
datrun.stimulus.y_end = 200;
datarun.field_height = 200;
datarun.field_width = 200;


%% plot individual RFs
cell_id = 7428;
cell_index = get_cell_indices(datarun,cell_id);
window_size = 30;
figure(2); clf;
image_scale_factor = [8,8];

image(matrix_scaled_up(norm_image(datarun.stas.rfs{cell_index}), image_scale_factor(1)))

%plot_rf(datarun,cell_id, 'scale', image_scale_factor(1))


[weights, selection, extras] = select_cone_weights(datarun, cell_id,...
                                            'thresh', 0.2,...
                                            'radius', [0 inf], 'polarity', 1,...
                                            'contiguity', true,'scale', 3.0);   

connected_cones = find(selection);

% color connected cones blue, others dark gray
cone_colors = 0.15 * ones(length(datarun.cones.centers(:,1)), 3);
cone_colors(connected_cones,:) = repmat([1 0 0], length(connected_cones),1);

% plot original cone mosaic
plot_cone_mosaic(datarun, 'fig_or_axes', 2, 'bg_color', [0.5 0.5 0.5], 'clear', false,...
                'cone_size', 15, 'cone_colors', cone_colors, 'scale', image_scale_factor,...
                'cone_roi', connected_cones);

rgc_center = datarun.stas.rf_coms{cell_index};
begin_x = rgc_center(1) - window_size;
begin_y = rgc_center(2) - 15;
end_x = rgc_center(1) + window_size;
end_y = rgc_center(2) + 45;

axis([begin_x end_x begin_y end_y]*image_scale_factor(1))
print(2, '~/Desktop/fig2.pdf', '-dpdf')

%%
cell_id = 31;

cell_index = get_cell_indices(datarun, cell_id);

figure(1); clf;
image(matrix_scaled_up(norm_image(datarun.stas.rfs{cell_index}), image_scale_factor(1)))

[weights, selection, extras] = select_cone_weights(datarun, cell_id,...
                                            'thresh', 0.2,...
                                            'radius', [0 inf], 'polarity', 1,...
                                            'contiguity', true,'scale', 3.0);   

connected_cones = find(selection);

% color connected cones blue, others dark gray
cone_colors = 0.15 * ones(length(datarun.cones.centers(:,1)), 3);
cone_colors(connected_cones,:) = repmat([0 1 0], length(connected_cones),1);

% plot original cone mosaic
plot_cone_mosaic(datarun, 'fig_or_axes', 1, 'bg_color', [0.5 0.5 0.5], 'clear', false,...
                'cone_size', 15, 'cone_colors', cone_colors, 'scale', image_scale_factor,...
                'cone_roi', connected_cones);
                
                
axis([begin_x end_x begin_y end_y]*image_scale_factor(1))
print(1, '~/Desktop/fig1.pdf', '-dpdf')



%%
cell_id = 7427;
cell_index = get_cell_indices(datarun,cell_id);

figure(3); clf;
image(matrix_scaled_up(norm_image(datarun.stas.rfs{cell_index}), image_scale_factor(1)))


[weights, selection, extras] = select_cone_weights(datarun, cell_id,...
                                            'thresh', 0.2,...
                                            'radius', [0 inf], 'polarity', 1,...
                                            'contiguity', true,'scale', 3.0);   

connected_cones = find(selection);

% color connected cones blue, others dark gray
cone_colors = 0.15 * ones(length(datarun.cones.centers(:,1)), 3);
cone_colors(connected_cones,:) = repmat([0 0 1], length(connected_cones),1);

plot_cone_mosaic(datarun, 'fig_or_axes', 3, 'bg_color', [0.5 0.5 0.5], 'clear', false,...
                'cone_size', 15, 'cone_colors', cone_colors, 'scale', image_scale_factor,...
                'cone_roi', connected_cones);

axis([begin_x end_x begin_y end_y]*image_scale_factor(1))
print(3, '~/Desktop/fig3.pdf', '-dpdf')








