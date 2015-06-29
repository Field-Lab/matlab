% runDS
% Colleen Rhoades
% June 2015
% rhoades@stanford.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The purpose of this file is to analyze the results from drifting gratings.

% You will see a graph showing the number of total spikes during the stimulus presentation versus the direction of the stimulus.
% If a direction is preferred by the cell, then the spike count should be higher.
% The plot is divided for different spatial and temporal frequencies.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear

%%%%%%%%%%%%%%%%%%%%%%%%%% INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load the Data for the DS run
run_opt.data_set = '2015-05-27-5';
run_opt.data_run = 'data003';

interval = 5; % each stimulus was displayed for 5 seconds

% Where to save the data
filepath= ['/Users/colleen/Desktop/DS/', run_opt.data_set, '/', run_opt.data_run, '/'];

% You can give the cells as all of them (datarun{2}.cell_ids) or give
% specific vision ids
% Find the cell to run by mapping a large cell EI from a white noise run
% cells = datarun{2}.cell_ids;
cells = 18;


%%%%%%%%%%%%%%%%%%%%%%%%%% END INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


datarun{1}.names.rrs_params_path=['/Volumes/Analysis/2015-05-27-5/', run_opt.data_run, '/', run_opt.data_run, '.params'];
datarun{2}.names.rrs_neurons_path=['/Volumes/Analysis/2015-05-27-5/', run_opt.data_run, '/', run_opt.data_run, '.neurons'];
datarun{2}.names.stimulus_path=['/Volumes/Data/2015-05-27-5/Visual/s', run_opt.data_run(end-1:end)];
opt=struct('verbose',1,'load_params',1,'load_neurons',1,'load_obvius_sta_fits',true);

datarun=load_data(datarun,opt);
datarun=map_cell_types(datarun, struct('map',[1 2],'verbose',true));

% If the trigger interval in LabView is set long enought (~6 seconds for 5
% second stimuli), then this trigger_iti_thr should be fine.
datarun{2}=load_stim(datarun{2},'correction_incomplet_run', 0, 'trigger_iti_thr', 0.1); % manually set threshold until got right number of triggers

spatial_tested = fliplr(datarun{2}.stimulus.params.SPATIAL_PERIOD);
temporal_tested = fliplr(datarun{2}.stimulus.params.TEMPORAL_PERIOD);


% get spikes per bin

% These triggers depend on the Labview trigger. 7.5 seconds prevents
% triggers from being missed on the large spatial scales
tr=datarun{2}.stimulus.triggers; % all start triggers

for x = 1:length(cells)
    x % see progress
    cell_to_run = cells(x);
    cell_indices1=get_cell_indices(datarun{1}, [cell_to_run]);
    cell_indices2=get_cell_indices(datarun{2},[cell_to_run]);
    
    % Get spikes times for the specified cell   
    n = datarun{2}.spikes{cell_indices2(1)}';
    psth_r = [];
    for z = 1:size(tr,2)
        ind = find(tr(z) == datarun{2}.triggers);
        spatial = datarun{2}.stimulus.trials(z).SPATIAL_PERIOD;
        temporal = datarun{2}.stimulus.trials(z).TEMPORAL_PERIOD;
        direction = datarun{2}.stimulus.trials(z).DIRECTION;
        num_triggers = interval/(temporal/120); % Number of times a band passed over one spot
        sub_spacing= datarun{2}.triggers(ind:ind+num_triggers) - datarun{2}.triggers(ind); % How to divide by the trial into the timing of the bands
        start = 0;
        stop = temporal/120;
        
        
        
        % first trial at the top
        for i=1:length(sub_spacing)-1,
            h=n-tr(z)-sub_spacing(i); % align the spike times
            stop = sub_spacing(i+1) - sub_spacing(i);
            hh=find(h>=start & h<=stop); % Find spikes for each grating that passed over the spot
            psth_r=[psth_r; (h(hh)*1000)', repmat(length(sub_spacing)-i,[length(hh),1]);];
            
            
        end
        psth_mat{z} = psth_r; % Put all the psths from all the temporal, spatial and direction combinations into a cell array
        psth_r = [];
        
    end
    
    % Organize the data into [direction, spatial, temporal, number of spikes in
    % 5 second trial]
    by_trial = zeros(size(tr,2),4);
    for i = 1:size(tr,2)
        by_trial(i,1) = datarun{2}.stimulus.trials(i).DIRECTION;
        by_trial(i,2) = datarun{2}.stimulus.trials(i).SPATIAL_PERIOD;
        by_trial(i,3) = datarun{2}.stimulus.trials(i).TEMPORAL_PERIOD;
        by_trial(i,4) = size(psth_mat{i},1); % number of spikes in that 5 seconds regardless of timing
        
    end
    
    % Plot the data separated out by either temporal or spatial period
    spat64 = find(by_trial(:,3) == temporal_tested(1));
    fig = figure;
    set(fig,'PaperPositionMode','auto')
    set(fig, 'PaperSize', [12 8], 'PaperPosition', [0 0 12,8]);
    set(fig, 'Position', [0 0 1200 800]);
    set(fig, 'PaperUnits', 'inches')
    set(fig, 'PaperPosition', [0 0 12 8]);
    
    
    for i = spat64
        subplot(2,3,1)
        plot(by_trial(i,1), by_trial(i,4), 'k.')
        hold on
    end
    title(['Temporal = ', num2str(temporal_tested(1))])
    xlabel('Direction')
    ylabel('Spikes count during stimulus')
    
    clear spat64
    
    spat64 = find(by_trial(:,3) == temporal_tested(2));
    for i = spat64
        subplot(2,3,2)
        
        plot(by_trial(i,1), by_trial(i,4), 'k.')
        hold on
    end
    
    title(['Temporal = ', num2str(temporal_tested(2))])
    xlabel('Direction')
    ylabel('Spikes count during stimulus')
    clear spat64
    
    
    spat64 = find(by_trial(:,3) == temporal_tested(3));
    for i = spat64
        subplot(2,3,3)
        
        plot(by_trial(i,1), by_trial(i,4), 'k.')
        hold on
    end
    
    title(['Temporal = ', num2str(temporal_tested(3))])
    xlabel('Direction')
    ylabel('Spikes count during stimulus')
    
    clear spat64
    
    spat64 = find(by_trial(:,2) == spatial_tested(1));
    for i = spat64
        subplot(2,3,4)
        
        plot(by_trial(i,1), by_trial(i,4), 'k.')
        hold on
    end
    title(['Spatial = ', num2str(spatial_tested(1))])
    xlabel('Direction')
    ylabel('Spikes count during stimulus')
    
    clear spat64
    
    spat64 = find(by_trial(:,2) == spatial_tested(2));
    
    for i = spat64
        subplot(2,3,5)
        
        plot(by_trial(i,1), by_trial(i,4), 'k.')
        hold on
    end
    title(['Spatial = ', num2str(spatial_tested(2))])
    xlabel('Direction')
    ylabel('Spikes count during stimulus')
    
    clear spat64
    
    
    spat64 = find(by_trial(:,2) == spatial_tested(3));
    
    for i = spat64
        subplot(2,3,6)
        
        plot(by_trial(i,1), by_trial(i,4), 'k.')
        hold on
    end
    title(['Spatial = ', num2str(spatial_tested(3))])
    xlabel('Direction')
    ylabel('Spikes count during stimulus')
    
    suptitle({run_opt.data_set;['Cell ', num2str(cell_to_run)]});
    
    
    if ~exist(filepath)
        mkdir(filepath)
    end
    
    % set(gcf, 'renderer', 'opengl')
    
    export_fig([filepath, 'Cell_',num2str(cell_to_run)], '-pdf')
    % print(gcf,[filepath, 'Cell_',num2str(cell_to_run)], '-depsc', '-opengl');
    close(gcf);
end
