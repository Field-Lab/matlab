path_and_name{1,1} = '/snle/lab/Experiments/Array/Analysis/2008-08-27-5/data003/data003/data003';
path_and_name{1,2} = 'plantain';
path_and_name{2,1} = '/snle/lab/Experiments/Array/Analysis/2008-04-22-5/data006/data006';
path_and_name{2,2} = 'plum';
path_and_name{3,1} = '/snle/lab/Experiments/Array/Analysis/2008-08-26-2/data001-s6369-s9551/data001-s6369-s9551';
path_and_name{3,2} = 'blueberry';
path_and_name{4,1} = '/snle/lab/Experiments/Array/Analysis/2008-12-12-1/data005/data005';
path_and_name{4,2} = 'butterfly';
path_and_name{5,1} = '/snle/lab/Experiments/Array/Analysis/2007-08-21-1/data003/data003';
path_and_name{5,2} = 'pomegranate';
path_and_name{6,1} = '/snle/lab/Experiments/Array/Analysis/2009-02-28-0/data006/data006';
path_and_name{6,2} = 'cherry';
path_and_name{7,1} = '/snle/lab/Experiments/Array/Analysis/2008-03-25-3/data002/data002';
path_and_name{7,2} = 'cherimoya';
path_and_name{8,1} = '/snle/lab/Experiments/Array/Analysis/2008-05-13-3/data006/data006';
path_and_name{8,2} = 'kiwi';
path_and_name{9,1} = '/snle/lab/Experiments/Array/Analysis/2008-04-30-2/data004/data004/data004';
path_and_name{9,2} = 'mango';
path_and_name{10,1} = '/snle/lab/Experiments/Array/Analysis/2007-03-27-2/data014/data014/data014';
path_and_name{10,2} = 'grapes';
path_and_name{10,1} = '/snle/lab/Experiments/Array/Analysis/2008-08-27-0/data001-s5661-s9260/data001-s5661-s9260';
path_and_name{10,2} = 'peach';



num_datasets = length(path_and_name(:,1));
verbose = 0;

mosaic_counter = 0;
clear mosaic_unique_cone_fraction mosaic_unique_cone_fraction_sd data_SD_PIs

mosaic_info = struct;
for dataset = 1:num_datasets
    m




rand('twister', 11111);
for dataset = 1:num_datasets
    clear datarun new_datarun sim_datarun
    datarun = load_data(path_and_name{dataset,1});
    datarun = load_params(datarun,struct('verbose',1));  
    datarun = load_sta(datarun,'load_sta',[]);
    datarun = import_single_cone_data(datarun, path_and_name{dataset, 2});
    
    num_on_midgets = length(datarun.cell_types{3}.cell_ids);
    num_off_midgets = length(datarun.cell_types{4}.cell_ids);
    
    if num_on_midgets > 20
        on_midget_flag = true;
    else
        on_midget_flag = false;
    end
    if num_off_midgets > 20
        off_midget_flag = true;
    else
        off_midget_flag = false;
    end
    
    if on_midget_flag && off_midget_flag
        cell_types = [3,4];
    elseif on_midget_flag && ~off_midget_flag
        cell_types = 3;
    elseif ~on_midget_flag && off_midget_flag
        cell_type = 4;
    else
        cell_types = [];
    end

    mosaic_info = setfield(mosaic_info, path_and_name{dataset, 1});
    
    for tp = 1:length(cell_types)
        cell_type = cell_types(tp);

        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % compute connectivity matrix
        [mosaic_weights, selection] = select_cone_weights(datarun, {cell_type}, 'thresh', 0.05, 'radius', [0 6], 'polarity', 1,...
                                                            'contiguity', true, 'scale', 2.5);                                        
        connectivity_two = mosaic_weights .* selection;

        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % compute surround connectivity
        [mosaic_weights, selection] = select_cone_weights(datarun, {cell_type}, 'thresh', 0.05, 'radius', [2 8], 'polarity', -1,...
                                                        'contiguity', true, 'scale', 2.0);                                        
        surround_connectivity = mosaic_weights .* selection;


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % compare center versus surround connectivity
        surround_volumes = sum(surround_connectivity);
        center_volumes = sum(connectivity_two);
        % ensure the ratios are real numbers
        kosher_indices = find(center_volumes > 0);
        surround_center_ratios = abs(surround_volumes(kosher_indices) ./ center_volumes(kosher_indices));


        hist(surround_center_ratios, [0:0.01:1.0])
        
        
        mosaic_counter = mosaic_counter + 1;
        
        mosaic_ratios(mosaic_counter) = median(surround_center_ratios);
        
        mosaic_info = setfield(
        
    end
end


hist(mosaic_ratios)
   
        
