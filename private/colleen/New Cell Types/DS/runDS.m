% ds_cells
run_opt.data_set = '2008-04-30-1';
run_opt.data_run = 3; % 12-19 for 2007-03-27, 2-11 for 2007-08-24, 13-17 for 2005-04-26
if strcmp(run_opt.data_set, '2008-04-30-1')
    datarun{1}.names.rrs_params_path='/Volumes/Analysis/2008-04-30-1/data009/data009-map/data009-from-data005.params';
    datarun{2}.names.rrs_neurons_path=('/Volumes/Analysis/2008-04-30-1/data009/data009-map/data009-map.neurons');
    datarun{2}.names.stimulus_path=('/Volumes/Archive/2005-01-21-0/Visual/s03');
end
% opt=struct('verbose',1,'load_params',1,'load_neurons',1,'load_obvius_sta_fits',true);
% datarun=load_data(datarun,opt);
% datarun=map_cell_types(datarun, struct('map',[1 2],'verbose',true));
datarun{2}=load_stim(datarun{2},'correction_incomplet_run', 0);

% k=1; kmin=1; kmax=size(datarun{2}.stimulus.trials,2); hk=loop_slider(k,kmin,kmax);
% 
% while k
%     if ~ishandle(hk)
%         break % script breaks until figure is closed
%     end
%     for z = 1:size(datarun{2}.stimulus.trials,2)
%         
%         k=round(get(hk,'Value'));
%         spatial = datarun{2}.stimulus.trials(z).SPATIAL_PERIOD
%         temporal = datarun{2}.stimulus.trials(z).TEMPORAL_PERIOD
%         direction = datarun{2}.stimulus.trials(z).DIRECTION
%         
%         
%         
%         
%         %         run_opt.config_num = findStimType(datarun, run_opt.contrast, run_opt.direction);
        cell_indices1=get_cell_indices(datarun{1}, [4802]);
        cell_indices2=get_cell_indices(datarun{2},[4802]);
%         cell_x_pos = cellfun( @(X) X.mean(1), datarun{1}.vision.sta_fits); % x axis position of all STA cells
%         [~, cell_sort_idx] = sort(cell_x_pos(cell_indices1)); % indicies of how to sort
%         
%         %cell_indices sorted by their x coordinate of the RF from the STA
%         cell_indices1 = cell_indices1(cell_sort_idx); % cell_indices1 is now indexes in order from lowest to highest firing rate
%         cell_indices2 = cell_indices2(cell_sort_idx);
%         
%         % Find trial start and stop times
%         %             start = 0;
%         %             stop = mean(datarun{2}.stimulus.triggers(2:2:end) - datarun{2}.triggers(1:2:end));
%         tr=datarun{2}.stimulus.triggers; % all start triggers
%         run_opt.config_num = 1;
%         t=find(datarun{2}.stimulus.trial_list==run_opt.config_num); %find the times when all the stimulus type 2 starts
%         tr=tr(1);
%         tr= 0:temporal/120:5;
%         start = 0;
%         stop = temporal/120;
%         
%         %         if size(run_opt.cell_type,2) == 1
%         %         cell_indices1 = [cell_indices1{1}];
%         %         cell_indices2 = [cell_indices2{1}];
%         %         %     cell_x_pos = [cell_x_pos{1}];
%         %         [~, cell_sort_idx] = sort(cell_x_pos(cell_indices1)); % indicies of how to sort
%         
%         %     else
%         %
%         %         cell_indices1 = [cell_indices1{1}, cell_indices1{2}];
%         %         cell_indices2 = [cell_indices2{1}, cell_indices2{2}];
%         %         %     cell_x_pos = [cell_x_pos{1}, cell_x_pos{2}];
%         %         [~, cell_sort_idx] = sort(cell_x_pos(cell_indices1)); % indicies of how to sort
%         %     end
%         
%         %cell_indices sorted by their x coordinate of the RF from the STA
%         cell_indices1= cell_indices1(cell_sort_idx); % cell_indices1 is now indexes in order from lowest to highest firing rate
%         cell_indices2 = cell_indices2(cell_sort_idx);
%         
%         
%         
%         run_opt.raster = 1;
%         if run_opt.raster %raster
%             %         cell_indices1 = [cell_indices1{1}, cell_indices1{2}];
%             %         cell_indices2 = [cell_indices2{1}, cell_indices2{2}];
%             %         cell_x_pos = [cell_x_pos{1}, cell_x_pos{2}];
%             %         [~, cell_sort_idx] = sort(cell_x_pos(cell_indices1)); % indicies of how to sort
%             
%             %cell_indices sorted by their x coordinate of the RF from the STA
%             cell_indices1= cell_indices1(cell_sort_idx); % cell_indices1 is now indexes in order from lowest to highest firing rate
%             cell_indices2 = cell_indices2(cell_sort_idx);
%             
%             
%             % Takes in start and stop time (0-0.7274)
%             % Spikes of the cell with the lowest firing rate first
%             % start time of each stimulus type trigger
%             % Finds the spikes that happened on a cell from stimulus onset to end
%             % Plot those spike times on the x axis versus the trial number on the y
%             % axis
%             % If tracking motion, the cell should respond to the bar at the same
%             % time on every trial
%             
%             psth_r = psth_raster(start,stop,datarun{2}.spikes{cell_indices2(1)}',tr);
%             
%             % Title is the cell id according to vision and the mean firing rate
%             title(sprintf('Cell ID: %d   Direction:%d   Temporal: %d   Spatial: %d', datarun{2}.cell_ids(cell_indices2(1)), direction, temporal, spatial ))
%             xlabel('ms')
%             drawnow
%                         uiwait;
%         end
%     end
% end

% get spikes per bin

     tr=datarun{2}.triggers; % all start triggers
     tr=tr(1:10:end-1); % all start triggers
        run_opt.config_num = 1;
%         t=find(datarun{2}.stimulus.trial_list==run_opt.config_num); %find the times when all the stimulus type 2 starts
%         tr=tr(1);

        
n = datarun{2}.spikes{cell_indices2(1)}';
psth_r = [];
for z = 1:size(tr,1)
        spatial = datarun{2}.stimulus.trials(z).SPATIAL_PERIOD;
        temporal = datarun{2}.stimulus.trials(z).TEMPORAL_PERIOD;
        direction = datarun{2}.stimulus.trials(z).DIRECTION;
                sub_spacing= 0:temporal/120:4.99;
        start = 0;
        stop = temporal/120;

        

        
for i=1:length(sub_spacing),
        h=n-tr(z)-sub_spacing(i);  
        hh=find(h>=start & h<=stop);
        psth_r=[psth_r; (h(hh)*1000)', repmat(length(sub_spacing)-i,[length(hh),1]);];
        
        
end
psth_mat{z} = psth_r;
psth_r = [];

end


