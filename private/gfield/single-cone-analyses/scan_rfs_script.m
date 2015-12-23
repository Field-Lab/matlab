path_and_name{1,1} = '/snle/lab/Experiments/Array/Analysis/2008-08-27-5/data003/data003/data003';
path_and_name{1,2} = 'plantain';

path_and_name{1,1} = '/snle/lab/Experiments/Array/Analysis/2008-08-27-0/data001-s5661-s9260/data001-s5661-s9260';
path_and_name{1,2} = 'peach';

path_and_name{1,1} = '/snle/lab/Experiments/Array/Analysis/2009-04-13-5/data005/data005';
path_and_name{1,2} = 'apricot';



path_and_name{1,1} = '/snle/lab/Experiments/Array/Analysis/2008-04-22-5/data006/data006';
path_and_name{1,2} = 'plum';

path_and_name{1,1} = '/snle/lab/Experiments/Array/Analysis/2008-05-13-3/data006/data006';
path_and_name{1,2} = 'kiwi';

path_and_name{1,1} = '/snle/lab/Experiments/Array/Analysis/2008-08-26-2/data001-s6369-s9551/data001-s6369-s9551';
path_and_name{1,2} = 'blueberry';

path_and_name{1,1} = '/Analysis/gfield/2008-04-22-5/data006-1800s-3600s/data006/data006';
path_and_name{1,2} = 'plum';


cell_type = 2;


% load data
datarun = load_data(path_and_name{1,1});
datarun = load_params(datarun,struct('verbose',1));  
datarun = load_sta(datarun,'load_sta',[]);
%datarun = import_single_cone_data(datarun, path_and_name{1, 2});


datarun = get_sta_summaries(datarun, {cell_type}, 'keep_stas', false);
datarun = set_polarities(datarun);

% portraits
plot_rf_portraits(datarun, {cell_type}, 'plot_radius', 25);





% get cell indicies
cell_indices = get_cell_indices(datarun, {cell_type});
num_cells = length(cell_indices);

% get and store significant stixels
sig_stix_params = struct('select', 'max', 'thresh', 2.0);
datarun = get_significant_stixels(datarun, {cell_type}, 'thresh_params', sig_stix_params);


% compute profiles
datarun = get_rfs(datarun, {cell_type},'sig_stixels', datarun.stas.significant_stixels);

plot_rf(datarun, datarun.cell_ids(cell_indices(1)))

% compute COMs
com_params = struct('exclude_outliers', true, 'distance_thresh', 50, 'positive', true);
datarun = get_rf_coms(datarun, {cell_type}, 'com_params', com_params, 'sig_stixels', datarun.stas.significant_stixels);






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute connectivity matrix
[mosaic_weights, selection] = select_cone_weights(datarun, {cell_type});                                        
connectivity = mosaic_weights .* selection;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute connectivity matrix
%[connectivity, new_datarun, connectivity_extras] = extract_connectivity(datarun, cell_type, 'remove_cones', 'SU', 'required_sign', 'positive',...
%                                            'min_convergence', RGC_min_convergence, 'max_convergence', RGC_max_convergence,...
%                                            'max_radius', radius_scaler, 'normalize', true);

[mosaic_weights, selection, extras] = select_cone_weights(datarun, {cell_type}, 'thresh', 0.05, 'radius', [0 inf], 'polarity', 1,...
                                            'contiguity', true, 'scale', 3.0, 'remove_cones', 'SU');                                        
connectivity_two = mosaic_weights .* selection;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute surround connectivity
[mosaic_weights, selection] = select_cone_weights(datarun, {cell_type}, 'thresh', 0.05, 'radius', [2 8], 'polarity', -1,...
                                                    'contiguity', true, 'scale', 3.0);                                        
surround_connectivity = mosaic_weights .* selection;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compare center versus surround connectivity
surround_volumes = sum(surround_connectivity);
center_volumes = sum(connectivity_two);
% ensure the ratios are real numbers
kosher_indices = find(center_volumes > 0);
surround_center_ratios = abs(surround_volumes(kosher_indices) ./ center_volumes(kosher_indices));

median(surround_center_ratios)
hist(surround_center_ratios, [0:0.02:0.5])



window_scale = 20;
for RGC = 1:length(cell_indices)
    
    % method 1
    figure(1)
    clf
    plot_rf(datarun, datarun.cell_ids(cell_indices(RGC)), 'foa', 1)
    hold on
    temp_rf_fit = datarun.cones.rf_fits{cell_indices(RGC)};
    rf_center = temp_rf_fit.center;
    
    % plot all cone locations
    for cone = 1:length(connectivity(:,RGC));
        plot(datarun.cones.centers(cone,1), datarun.cones.centers(cone, 2), 'wo', 'MarkerFaceColor', 'w')
    end
    
    % plot connected cone locations
    temp_cone_indices = find(connectivity(:, RGC) > 0);
    for cone = 1:length(temp_cone_indices)
        plot(datarun.cones.centers(temp_cone_indices(cone),1), datarun.cones.centers(temp_cone_indices(cone), 2), 'ro', 'MarkerFaceColor', 'r')
    end

    plot(rf_center(1), rf_center(2), 'ko', 'MarkerFaceColor', 'k')
    axis_window = [(rf_center(1)-window_scale) (rf_center(1)+window_scale) (rf_center(2)-window_scale) (rf_center(2)+window_scale)];
    axis(axis_window)    
    hold off

    
    % method 2
    figure(2)
    clf
    plot_rf(datarun, new_datarun.cell_ids(cell_indices(RGC)), 'foa', 2)
    hold on
    temp_rf_fit = new_datarun.cones.rf_fits{cell_indices(RGC)};
    rf_center = temp_rf_fit.center;
    
    % plot all cone locations
    for cone = 1:length(connectivity_two(:,RGC));
        plot(new_datarun.cones.centers(cone,1), new_datarun.cones.centers(cone, 2), 'wo', 'MarkerFaceColor', 'w')
    end
    
    % plot connected cone locations
    temp_cone_indices = find(connectivity_two(:, RGC) > 0);
    for cone = 1:length(temp_cone_indices)
        plot(new_datarun.cones.centers(temp_cone_indices(cone),1), new_datarun.cones.centers(temp_cone_indices(cone), 2), 'ro', 'MarkerFaceColor', 'r')
    end

    plot(rf_center(1), rf_center(2), 'ko', 'MarkerFaceColor', 'k')
    axis_window = [(rf_center(1)-window_scale) (rf_center(1)+window_scale) (rf_center(2)-window_scale) (rf_center(2)+window_scale)];
    axis(axis_window)    
    hold off    
    
    drawnow
    pause
    
end

    
                                        
                                          

