run_opt.cell_types = {'Off midget', 'Off parasol', 'On midget', 'On parasol'};

run_opt.cell_type = 'On parasol';

datarun{1}.names.rrs_params_path='/Volumes/Analysis/2007-01-23-5/data010/data010.params';
datarun{2}.names.rrs_neurons_path='/Volumes/Analysis/2007-01-23-5/data010/data010.neurons';
%             datarun{2}.names.stimulus_path=sprintf('/Volumes/Analysis/2007-03-27-1/stimuli/s%02d', run_opt.data_run);

opt=struct('verbose',1,'load_params',1,'load_neurons',1,'load_obvius_sta_fits',true, 'load_all', true);
datarun=load_data(datarun,opt);
datarun=map_cell_types(datarun, struct('map',[1 2],'verbose',true));
%         datarun{2}=load_stim(datarun{2},'correction_incomplet_run', 0);


% Get the right configuration number
%     run_opt.config_num = findStimType(datarun, run_opt.contrast, run_opt.direction);

% Get indicies for each cell type

    % Gets the indicies used by vision of the particular cell type
        
        
        % Get indices for specified cell type and order by RF position
     cell_indices1=get_cell_indices(datarun{1},{run_opt.cell_type});
    cell_indices2=get_cell_indices(datarun{2},{run_opt.cell_type});
    
    cell_x_pos = cellfun( @(X) X.mean(1), datarun{1}.vision.sta_fits); % x axis position of all STA cells
    [~, cell_sort_idx] = sort(cell_x_pos(cell_indices1)); % x axis position of only on midget cells, indexes of how to sort
    
    %cell_indices sorted by their x coordinate of the RF from the STA
    cell_indices1 = cell_indices1(cell_sort_idx); % cell_indices1 is now indexes in order from lowest to highest firing rate
    cell_indices2 = cell_indices2(cell_sort_idx);
    
    % Find trial start and stop times
    start = 0;
    stop = mean(datarun{2}.triggers(2:2:end) - datarun{2}.triggers(1:2:end));
    tr=datarun{2}.triggers(1:2:end); % all start triggers
    t=find(datarun{2}.stimulus.trial_list==run_opt.config_num); %find the times when all the stimulus type 2 starts
    tr=tr(t);
    

% Grouped pooled together
if size(run_opt.cell_type,2) == 1
    cell_indices1 = [cell_indices1{1}];
    cell_indices2 = [cell_indices2{1}];
    %     cell_x_pos = [cell_x_pos{1}];
    [~, cell_sort_idx] = sort(cell_x_pos(cell_indices1)); % indicies of how to sort
    
else
    
    cell_indices1 = [cell_indices1{1}, cell_indices1{2}];
    cell_indices2 = [cell_indices2{1}, cell_indices2{2}];
    %     cell_x_pos = [cell_x_pos{1}, cell_x_pos{2}];
    [~, cell_sort_idx] = sort(cell_x_pos(cell_indices1)); % indicies of how to sort
end

%cell_indices sorted by their x coordinate of the RF from the STA
cell_indices1= cell_indices1(cell_sort_idx); % cell_indices1 is now indexes in order from lowest to highest firing rate
cell_indices2 = cell_indices2(cell_sort_idx);