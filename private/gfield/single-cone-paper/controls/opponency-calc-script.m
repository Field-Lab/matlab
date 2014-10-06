
% script for comparing the surround strength with the center strength and correlating this with purity

% get data names and paths
[LMS_paths, LMS_names, ] = get_LMS_paths('high');

num_datasets = length(LMS_paths);
verbose = 0;
mosaic_counter = 0;
cell_cutoff = 60;

for dataset = 1:num_datasets
    clear datarun 
    datarun = load_data(LMS_paths{dataset});
    datarun = load_params(datarun,struct('verbose',1));  
    datarun = load_sta(datarun,'load_sta',[]);
    datarun = import_single_cone_data(datarun, LMS_names{dataset});

    num_on_midgets = length(datarun.cell_types{3}.cell_ids);
    num_off_midgets = length(datarun.cell_types{4}.cell_ids);

    if num_on_midgets > cell_cutoff
        on_midget_flag = true;
    else
        on_midget_flag = false;
    end
    if num_off_midgets > cell_cutoff
        off_midget_flag = true;
    else
        off_midget_flag = false;
    end
    
    if on_midget_flag && off_midget_flag
        cell_types = [3,4];
    elseif on_midget_flag && ~off_midget_flag
        cell_types = 3;
    elseif ~on_midget_flag && off_midget_flag
        cell_types = 4;
    else
        cell_types = [];
    end

    for tp = 1:length(cell_types)
        mosaic_counter = mosaic_counter + 1;
        cell_type = cell_types(tp);

       % ss_params.thresh = 2.5;
       % distance_threshold = 40;
        temp_cell_indices = get_cell_indices(datarun, {cell_type});
       % for RGC = 1:length(temp_cell_indices);
       %     RGC
       %     temp_sta = get_sta(datarun, datarun.cell_ids(temp_cell_indices(RGC)));
       %      datarun.stas.rf_coms{temp_cell_indices(RGC)} = rf_com(temp_sta,'ss_params',ss_params, 'distance_thresh', distance_threshold, 'positive', true);
       % end
        % new fits
        datarun = fit_cone_rfs(datarun, {cell_type},'verbose', true, 'show_fit_params', false, 'foa_profile', 1, 'foa_2d', 2);
        temp_filename = [LMS_names{dataset},'-',num2str(tp)];
        save(temp_filename, 'datarun')

          %temp_filename = [LMS_names{dataset},'-',num2str(tp)];
          %load(temp_filename)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % extract center connectivity
        [center_weights, center_selection, extras] = select_cone_weights(datarun, {cell_type}, 'thresh', 2.5,...
                                                    'radius', [0 6], 'polarity', 0,...
                                                    'contiguity', false, 'scale', 3.0, 'remove_cones', 'SU');   

        center_connectivity = center_weights .* center_selection;
        new_datarun = extras.new_datarun;
        [num_cones, num_RGCs] = size(center_connectivity);                                     

        % compute purity of cone mosaic
        purity_indices = compute_opponency_index(center_connectivity, new_datarun.cones.types);
        purity_sd = std(purity_indices);
        mosaic_purities(mosaic_counter) = purity_sd;

    end
end
