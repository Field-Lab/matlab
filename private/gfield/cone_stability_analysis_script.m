%% define data to analyze 

% 2011-10-25-9
datarunA = load_data('/snle/lab/Experiments/Array/Analysis/2011-10-25-9/streamed/data006-0/data006-0');
datarunB = load_data('/snle/lab/Experiments/Array/Analysis/2011-10-25-9/streamed/data010-0/data010-0');
cd /snle/acquisition/maps/2011-10-25-9_f06_vorcones
load('map-0000.txt')
cone_map = map_0000; clear map_0000;


% 2011-12-13-2 d14-d18
datarunA = load_data('/snle/lab/Experiments/Array/Analysis/2011-12-13-2/streamed/data014-0/data014-0');
datarunB = load_data('/snle/lab/Experiments/Array/Analysis/2011-12-13-2/streamed/data018-0/data018-0');
cd /snle/acquisition/maps/2011-12-13-2_f14_vorcones
load('map-0000.txt')
cone_map = map_0000; clear map_0000;

% 2011-12-13-2 d18-d23
datarunA = load_data('/snle/lab/Experiments/Array/Analysis/2011-12-13-2/streamed/data018-0/data018-0');
datarunB = load_data('/snle/lab/Experiments/Array/Analysis/2011-12-13-2/streamed/data023-0/data023-0');
cd /snle/acquisition/maps/2011-12-13-2_f18_vorcones
load('map-0000.txt')
cone_map = map_0000; clear map_0000;



%% load data
datarunA = load_sta(datarunA, 'load_sta', []);
datarunA = load_params(datarunA);
datarunA = load_cones(datarunA);
datarunA.cones.mosaic = make_mosaic_struct(datarunA.cones.centers);

datarunB = load_sta(datarunB, 'load_sta', []);
datarunB = load_params(datarunB);
datarunB = load_cones(datarunB);
datarunB.cones.mosaic = make_mosaic_struct(datarunB.cones.centers);


cell_types = {1,2,3,4,5};
datarunA = get_sta_summaries(datarunA, 'all');
datarunB = get_sta_summaries(datarunB, 'all');
datarunA = get_sta_fits_from_vision(datarunA, 'all');
datarunB = get_sta_fits_from_vision(datarunB, 'all');
% map RGCs
temp_cell_type = {2};
cell_indices = get_cell_indices(datarunA, temp_cell_type);

for cn = length(cell_indices):-1:1
    figure(cn)
    plot_rf_summaries(datarunA, datarunA.cell_types{temp_cell_type{:}}.cell_ids(cn),...
                'foa',cn,'clear',true,'label',true,'plot_fits',true,'fit_color','k')
            
    plot_rf_summaries(datarunB, temp_cell_type,...
                'foa',cn,'clear',false,'label',true, 'label_color','r','plot_fits',true,'fit_color','r')
end


%% 2011-10-25-9 d06-d10
%off midgets
mapped_ids = [7353 213 413 424 438 511 572 753 756 785 1159 961 1068 932 1411 1431 1639 1817,...
                1876 2000 2074 -1 1909 2191 2438 2971 2881 3181 3076 3408 3683 3619 3859 4038,...
                3602 4068 3998 4055 4173 4428 4339 4486 4487 4682 4670 4833 4596 4806 -1 4967,...
                5120 5255 5011 5358 5163 5585 5481 -1 -1 -1 5765 6276 6098 6290 6381 6485 6706,...
                -1 6814 6890 6781 7113 7321 7381 -1 7666 7730 7759];
% off parasol
mapped_ids = [781 301 4774 1711 2086 3842 -1 3541 -1 4069 4666 4966 5327 5866 6091 6228,...
                5761 6812 7246 7471 7726];
% on midget
mapped_ids = [-1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 2194 2312 -1 2566 -1 -1 -1 -1,...
                3213 3291 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1];
            

%% 2011-12-13-2 d14-d18        
            
% on midget
mapped_ids = [-1 601 -1 859 -1 1518 1446 -1 -1 -1 -1 -1 -1 -1 -1 2269 2300 2362 2661,...
               2616 -1 -1 2821 2822 3035 -1 3083 3122 3245 3200 -1 -1 3648 3817 3857 3965 -1,...
               -1 -1 4726 4883 5300 5476 5719 -1 -1 7340 7471 -1];

% off midget
mapped_ids = [-1 602 766 811 916 753 1025 1024 -1 1321 1339 1351 -1 1564 1591 1743 1906,...
                -1 2091 2088 2251 -1 -1 -1 2581 -1 -1 2658 2686 2836 2929 2971 3031 -1,...
                3181 -1 -1 3244 -1 3331 2836 3392 3406 3541 3586 3632 -1 -1 -1 -1 3918 -1,...
                -1 4141 4171 4068 4216 4472 -1 4576 4577 3016 4861 5161 -1 -1 5357 -1 7036];

% off parasol
mapped_ids = [1126 2313 2881 3286 3361 3631 3901 3991 4291 4966 4996 6603 -1];
            
%% get voronoi info
datarunA.cones.mosaic.voronoi_masks_scale = 3; % For mapping based on 2x2 stixel run
datarunA = make_mosaic_struct(datarunA);
datarunA = make_voronoi_masks(datarunA);

%% Get cone size

% get cone locations
temp_locations = datarunA.cones.centers;
temp_weights = datarunA.cones.weights;
% get the index to the RGC that samples each cone strongest
[max_weights, max_weight_indices] =  max(temp_weights, [], 2);

window_size = 15;

summed_rf = zeros((window_size*2) +3);
for cn = 1:size(temp_locations,1);
    temp_rf = matrix_scaled_up(datarunA.stas.rfs{max_weight_indices(cn)}, 3);
    
    %imagesc(norm_image(temp_rf))
    
    x_begin = round(temp_locations(cn,1)*3) - (window_size + 1);
    y_begin = round(temp_locations(cn,2)*3) - (window_size + 1);
    x_end = round(temp_locations(cn,1)*3) + (window_size + 1);
    y_end = round(temp_locations(cn,2)*3) + (window_size + 1);

    cut_rf = temp_rf(y_begin:y_end, x_begin:x_end);
    %imagesc(norm_image(cut_rf))
    %pause
    summed_rf = cut_rf + summed_rf;
end

%summed_rf = summed_rf ./ size(temp_locations,1);

imagesc(norm_image(summed_rf))
    
plot(summed_rf(16,:))

% subtract baseline
rf_profile = summed_rf(16,:);
rf_profile = rf_profile - rf_profile(10);

x_points = [-5:1:5];

coef = [3,1];
fit_coef = nlinfit([-5:1:5], rf_profile(11:21), 'fit_Gaussian', coef)

plot(x_points, fit_Gaussian(fit_coef, x_points), 'r', [-15:1:17], rf_profile(:), 'k')

cone_size = fit_coef(2) * 2; % 2-sigma boundary;


%% loop over cells of a choosen type:

cell_type = {2};
cell_indices = get_cell_indices(datarunA, cell_type);
image_scale_factor = [600./datarunA.stimulus.field_height,600./datarunA.stimulus.field_height];
%image_scale_factor = [1 1];

% size of frame for plots (smaller for midgets, bigger for parasols)
window_size = 25 * image_scale_factor(1);

cd ~/Desktop/stability/

mask_size = 5;

for cc = 1:length(cell_indices)
    
    if mapped_ids(cc) == -1
        disp('cell did not map \n')
        continue
    else

        % get center of RGC RF
        rgc_center = datarunA.stas.fits{cell_indices(cc)}.mean;
        x_begin = (rgc_center(1) * image_scale_factor(1)) - window_size;
        y_begin = (rgc_center(2) * image_scale_factor(1)) - window_size;
        x_end = (rgc_center(1) * image_scale_factor(1)) + window_size;
        y_end = (rgc_center(2) * image_scale_factor(1)) + window_size;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        figure(1); clf
        %temp_pol = datarunA.stas.polarities{cell_indices(cc)};
        temp_pol = -1;
        image(norm_image(matrix_scaled_up(temp_pol*datarunA.stas.rfs{cell_indices(cc)},image_scale_factor(1))))
        axis square
        hold on

        plot_voronoi(datarunA, 'mode', 'manhattan', 'mask_space_self_max', mask_size, 'scale_plot', image_scale_factor(1))


        % extract connectivity
        [mosaic_weights, selection, extras] = select_cone_weights(datarunA, datarunA.cell_ids(cell_indices(cc)),...
                                                    'thresh', 0.1,...
                                                    'radius', [0 inf], 'polarity', 1,...
                                                    'contiguity', true,'scale', 3.0);   
        connected_cones = find(selection);

        % color connected cones blue, others dark gray
        cone_colors = 0.35 * ones(length(datarunA.cones.centers(:,1)), 3);
        cone_colors(connected_cones,:) = repmat([0 0 1], length(connected_cones),1);

        % plot original cone mosaic
        plot_cone_mosaic(datarunA, 'fig_or_axes', 1, 'bg_color', [], 'clear', false,...
                        'cone_size', 15, 'cone_colors', cone_colors, 'scale', image_scale_factor,...
                        'cone_circles', true, 'cone_radius', cone_size,...
                        'cone_roi', connected_cones);

        axis([x_begin x_end y_begin y_end])
        title(['cell id: ' num2str(datarunA.cell_types{cell_type{:}}.cell_ids(cc))])


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        figure(2); clf
        temp_index = get_cell_indices(datarunB, mapped_ids(cc));
        image(norm_image(matrix_scaled_up(temp_pol*datarunB.stas.rfs{temp_index},image_scale_factor(1))))
        axis square
        hold on
        plot_voronoi(datarunA, 'mode', 'manhattan', 'mask_space_self_max', mask_size, 'scale_plot', image_scale_factor(1))

        % plot later cone mosaic            
        plot_cone_mosaic(datarunB, 'fig_or_axes', 2, 'bg_color', [], 'clear', false,...
                        'cone_size', 6, 'cone_colors', [0.8 0.8 0.8], 'scale', image_scale_factor,...
                        'cone_circles', true, 'cone_radius', cone_size)    
        axis([x_begin x_end y_begin y_end])
        title(['cell id: ' num2str(mapped_ids(cc))])
        
        %plot voronoi tesselation
    %     for ptch = 1:num_patchs
    %         patch(test_stats(ptch).ConvexHull(:,1), test_stats(ptch).ConvexHull(:,2),'r','FaceAlpha', 0.1)
    %     end
    %     

        figure(3); clf; hold on
        image(ones(600,600,3) * 0.5);
        plot_voronoi(datarunA, 'mode', 'manhattan', 'mask_space_self_max', mask_size, 'scale_plot', image_scale_factor(1))
        plot_cone_mosaic(datarunA, 'fig_or_axes', 3, 'bg_color', [0.5 0.5 0.5], 'clear', false,...
                        'cone_size', 15, 'cone_colors', cone_colors, 'scale', image_scale_factor,...
                        'cone_circles', true, 'cone_radius', cone_size,...
                        'cone_roi', connected_cones);
            % plot later cone mosaic            
        plot_cone_mosaic(datarunB, 'fig_or_axes', 3, 'bg_color', [0.5 0.5 0.5], 'clear', false,...
                            'cone_size', 6, 'cone_colors', [0.8 0.8 0.8], 'scale', image_scale_factor,...
                            'cone_circles', true, 'cone_radius', cone_size)    
        axis([x_begin x_end y_begin y_end])
        axis square


        % plot RFs
    %     datarunA = get_sta_fits_from_vision(datarunA, datarunA.cell_ids(cell_indices(cc)));
    %     plot_rf_summaries(datarunA, datarunA.cell_ids(cell_indices(cc)), 'foa', 1, 'clear', false,...
    %                 'plot_fits', true, 'fit_color', [1 1 0], 'fit_width', 1,...
    %                 'label', true, 'scale', image_scale_factor(1),...
    %                 'label_color', [1 1 0], 'label_size', 16);
    %     

    %     plot_title = ['Cell ', num2str(datarunA.cell_ids(cell_indices(cc)))];
    %     title(plot_title)
          mkdir(num2str(datarunA.cell_ids(cell_indices(cc))))
          datarunA.cell_ids(cell_indices(cc))
          save_path = ['~/Desktop/stability/',num2str(datarunA.cell_ids(cell_indices(cc))),'/A'];
          print(1, save_path, '-dpdf')
          save_path = ['~/Desktop/stability/',num2str(datarunA.cell_ids(cell_indices(cc))),'/B'];
          print(2, save_path, '-dpdf')
          save_path = ['~/Desktop/stability/',num2str(datarunA.cell_ids(cell_indices(cc))),'/C'];
          print(3, save_path, '-dpdf')

        pause(0.5)
    end
end




%% plot RGC mosaic of choice over the top


plot_rf_summaries(datarunA, [3151 3138 2313], 'foa', 2, 'clear', true, ...
                'plot_fits', true, 'fit_color', [1 0 0], 'fit_width', 1,...
                'label', true, 'scale', image_scale_factor(1),...
                'label_color', [1 0 0], 'label_size', 14)


plot_rf_summaries(datarunA, [1336 1591 1741 2191 2251 2566 2656 2971 3031 3247 3331 3406 3586 4576 4621],...
                     'foa', 2, 'clear', false, ...
                     'plot_fits', true, 'fit_color', [0 1 0], 'fit_width', 1,...
                'label', true, 'scale', image_scale_factor(1),...
                'label_color', [0 1 0], 'label_size', 14)

plot_rf_summaries(datarunA, [601 858 2298 2821 3121 3245 3963], 'foa', 2, 'clear', false, ...
                'plot_fits', true, 'fit_color', [0 0 1], 'fit_width', 1,...
                'label', true, 'scale', image_scale_factor(1),...
                'label_color', [0 0 1], 'label_size', 14)
