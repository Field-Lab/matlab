

%%%% select cell

% OFF cells
%cell_id = 1099; % M off, L ON
cell_id = 3136; % L OFF, M, OFF
% On cells
%cell_id = 813; % L ON, M OFF

cell_id = 1381


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

cell_threshold = 0.03;

% plot background
% get size and color
y = datarun.stimulus.field_height;
x = datarun.stimulus.field_width;
rgb = [0.5 0.50 0.5];
% generate background matrix
plot_mat = cat(3,repmat(rgb(1),y,x),repmat(rgb(2),y,x),repmat(rgb(3),y,x));
% plot it
figure(10);clf;image(plot_mat);axis image; hold on


% plot_center and surround spiders
plot_cell_sampling(datarun,cell_id,'type','spider','fig_or_axes',10,'clear',0,'plot_cones',0,...
   'line_width',[realmin 2.5],'cell_colors',[0 0 0],'plot_radius',[0 6],'thresh',cell_threshold,'polarity',0,'contiguity',true)

type_colors = [1 0 0; 0 1 0; 0 0 1; 0 0 0];

% plot all cones
plot_cone_mosaic(datarun,'fig_or_axes',10,'cone_size',4,'clear',0,'bg_color',[], 'type_colors', type_colors)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get center information and plot center cones
% get center cones
clear selection_params
selection_params.radius = [0 6];
selection_params.contiguity = true;
selection_params.thresh = cell_threshold;
%selection_params.polarity = 1;
[weights, center_selection, extras] = select_cone_weights(datarun, cell_id, selection_params);

% plot center cone halos
cone_colors = [0.5 0.5 0.5] + (0.5 * repmat(datarun.stas.polarities{cell_index}, 1,3));
cone_colors = repmat(cone_colors,4,1);

plot_cone_mosaic(datarun,'fig_or_axes',10,'cone_size',7,'clear',0,'bg_color',[],...
                'cone_roi', center_selection, 'roi_highlight', false, 'type_colors', cone_colors)

% plot center cones
plot_cone_mosaic(datarun,'fig_or_axes',10,'cone_size',5,'clear',0,'bg_color',[],...
                'cone_roi', center_selection, 'roi_highlight', false, 'type_colors', type_colors)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get center information and plot center cones
% get center cones
clear selection_params
selection_params.radius = [0 6];
selection_params.contiguity = true;
selection_params.thresh = cell_threshold;
selection_params.polarity = -1;
[weights, center_selection, extras] = select_cone_weights(datarun, cell_id, selection_params);

% plot center cone halos
cone_colors = [0.5 0.5 0.5] + (-0.5 * repmat(datarun.stas.polarities{cell_index}, 1,3));
cone_colors = repmat(cone_colors,4,1);

plot_cone_mosaic(datarun,'fig_or_axes',10,'cone_size',7,'clear',0,'bg_color',[],...
                'cone_roi', center_selection, 'roi_highlight', false, 'type_colors', cone_colors)

% plot center cones
plot_cone_mosaic(datarun,'fig_or_axes',10,'cone_size',5,'clear',0,'bg_color',[], 'cone_roi', center_selection, 'roi_highlight', false)

print(10,'/snle/home/gfield/Desktop/spider','-dpdf')




