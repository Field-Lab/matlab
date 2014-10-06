%datarunA = load_data('/snle/acquisition/2010-09-24-1/data006/data006');
datarunA = load_data('/snle/lab/Experiments/Array/Analysis/2011-10-25-9/streamed/data006/data006');
datarunA = load_sta(datarunA, 'load_sta', []);
datarunA = load_params(datarunA);
datarunA = import_single_cone_data(datarunA, 'orange-6');

%datarunB = load_data('/snle/acquisition/2010-09-24-1/data035/data035');
datarunB = load_data('/snle/lab/Experiments/Array/Analysis/2011-10-25-9/streamed/data010/data010');
datarunB = load_sta(datarunB, 'load_sta', []);
datarunB = load_params(datarunB);
datarunB = import_single_cone_data(datarunB, 'orange-35');


plot_cone_mosaic(datarunB, 'fig_or_axes', 1,...
                'cone_size', 6, 'cone_colors', [0 0 0])
print(1, '~/Desktop/cone-mosaic-a.pdf','-dpdf')

plot_cone_mosaic(datarunA, 'fig_or_axes', 2,...
                'cone_size', 3, 'cone_colors', [1 1 1])
print(2, '~/Desktop/cone-mosaic-b.pdf','-dpdf')



%plot voronoi tesselation

test_stats = regionprops(cone_map,'centroid', 'convexhull');

num_patchs = length(test_stats);

figure(3); clf; hold on
for ptch = 1:num_patchs
    patch(test_stats(ptch).ConvexHull(:,1), test_stats(ptch).ConvexHull(:,2), 'r')
end
axis([1 320 1 320])
print(3,'~/Desktop/voronoi.pdf', '-dpdf')






datarun = load_data('2010-03-05-2', 'rf-13-apple-gf');
datarun = load_index(datarun);
datarun = load_sta(datarun, 'load_sta', []);
datarun = load_params(datarun);
datarun = import_single_cone_data(datarun);


plot_cone_mosaic(datarun, 'fig_or_axes', 1,...
                'cone_size', 3)

print(1,'~/Desktop/colored-cone.pdf', '-dpdf')


%%
datarunC = load_data('/snle/lab/Experiments/Array/Analysis/2010-09-24-1/data035/data035');
datarunC = load_sta(datarunC, 'load_sta', []);
datarunC = load_params(datarunC);

datarunC = get_sta_fits_from_vision(datarunC, 'all');
plot_rf_summaries(datarunC, {2},'foa', 3, 'plot_fits', true, 'label', true) 

print(3, '~/Desktop/off-parasol-rfs.pdf', '-dpdf')

%%
% load data
datarun = load_data(path_and_name{1,1});
datarun = load_params(datarun,struct('verbose',1));  
datarun = load_sta(datarun,'load_sta',[]);
datarun = set_polarities(datarun);
datarun = import_single_cone_data(datarun, path_and_name{1,2});    
datarun.cones.mosaic = make_mosaic_struct(datarunA.cones.centers);


cell_spec = {3};
temp_polarity = 1;

% plot background

% get size and color
y = datarunA.stimulus.field_height;
x = datarunA.stimulus.field_width;
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
[weights, center_selection, extras] = select_cone_weights(datarunA, cell_spec, selection_params);

summed_selection = sum(center_selection, 2);
sampled_indices = find(summed_selection > 0);
cone_center_roi = zeros(length(summed_selection),1);
cone_center_roi(sampled_indices) = 1;

% plot spiders
plot_cell_sampling(datarun,cell_spec,'type','spider','fig_or_axes',10,'clear',0,'plot_cones',0,...
   'line_width',[realmin 2.0],'cell_colors',...
    [0.5 0.5 0.5] + (0.5*temp_polarity),'plot_radius',[0 6],...
    'thresh', selection_params.thresh)

