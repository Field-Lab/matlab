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
run_opt.data_set = '2015-08-17-5';
run_opt.data_run = 'd01-29-norefit/data023';
location = 'Data';

interval = 5; % number of seconds each stimulus was displayed for 

% Where to save the data
filepath= ['/Users/colleen/Desktop/DS/', run_opt.data_set, '/', run_opt.data_run, '/'];

% You can give the cells as all of them (datarun{2}.cell_ids) or give
% specific vision ids
% Find the cell to run by mapping a large cell EI from a white noise run
cells = 'all'; % 'all' or give a vector of vision cell ids


%%%%%%%%%%%%%%%%%%%%%%%%%% END INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

datarun{1}.names.rrs_params_path=['/Volumes/Analysis/', run_opt.data_set, '/', run_opt.data_run, '/', 'data023', '.params'];
% datarun{1}.names.rrs_sta_path=['/Volumes/Analysis/', run_opt.data_set, '/', run_opt.data_run, '/', 'data022', '.sta'];

datarun{2}.names.rrs_neurons_path=['/Volumes/Analysis/', run_opt.data_set, '/', run_opt.data_run, '/', 'data023', '.neurons'];
datarun{2}.names.stimulus_path=['/Volumes/',location,'/', run_opt.data_set, '/Visual/s', run_opt.data_run(end-1:end)];


% datarun{1}.names.rrs_params_path=['/Volumes/Analysis/', run_opt.data_set, '/', run_opt.data_run, '/', run_opt.data_run, '.params'];
% datarun{1}.names.rrs_sta_path=['/Volumes/Analysis/', run_opt.data_set, '/', run_opt.data_run, '/', run_opt.data_run, '.sta'];
% 
% datarun{2}.names.rrs_neurons_path=['/Volumes/Analysis/', run_opt.data_set, '/', run_opt.data_run, '/', run_opt.data_run, '.neurons'];
% datarun{2}.names.stimulus_path=['/Volumes/',location,'/', run_opt.data_set, '/Visual/s', run_opt.data_run(end-1:end)];
opt=struct('verbose',1,'load_sta', 1,'load_params',1,'load_neurons',1,'load_obvius_sta_fits',true);

datarun=load_data(datarun,opt);
datarun=map_cell_types(datarun, struct('map',[1 2],'verbose',true));

% If the trigger interval in LabView is set long enought (~6 seconds for 5
% second stimuli), then this trigger_iti_thr should be fine.
datarun{2}=load_stim(datarun{2},'correction_incomplet_run', 0, 'trigger_iti_thr', 0.1); % manually set threshold until got right number of triggers

if strcmp(cells, 'all')
    cells = datarun{2}.cell_ids;
end


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
    for z = 1:size(tr,2) % skip the first 12 trials due to contrast adapation
        ind = find(tr(z) == datarun{2}.triggers);
        spatial = datarun{2}.stimulus.trials(z).SPATIAL_PERIOD;
        temporal = datarun{2}.stimulus.trials(z).TEMPORAL_PERIOD;
        direction = datarun{2}.stimulus.trials(z).DIRECTION;
        rgb = datarun{2}.stimulus.trials(z).RGB;
           
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
%                if rgb(1) == 0.48
%                    psth_mat{z} = [];
%                end
               
        psth_r = [];
        
    end
    

%     for i = 25:size(tr,2)
%         figure
%         plot(psth_mat{i}(:,1), psth_mat{i}(:,2), 'o')
%     end
    
    
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
    fig = figure;
    set(fig,'PaperPositionMode','auto')
    set(fig, 'PaperSize', [12 8], 'PaperPosition', [0 0 12,8]);
    set(fig, 'Position', [0 0 1200 800]);
    set(fig, 'PaperUnits', 'inches')
    set(fig, 'PaperPosition', [0 0 12 8]);
    
    
    for t = 1:length(temporal_tested)
            spat64 = find(by_trial(:,3) == temporal_tested(t));

        for i = spat64
            subplot(2,max(length(temporal_tested), length(spatial_tested)),t)
            plot(by_trial(i,1), by_trial(i,4), 'k.')
            hold on
        end
    title(['Temporal = ', num2str(temporal_tested(t))])
    xlabel('Direction')
    ylabel('Spikes count during stimulus')
    
    clear spat64
    
    end
    
   
   for t = 1:length(spatial_tested)

    spat64 = find(by_trial(:,2) == spatial_tested(t));
    for i = spat64
        subplot(2,max(length(temporal_tested), length(spatial_tested)),t + max(length(temporal_tested), length(spatial_tested)))
        
        plot(by_trial(i,1), by_trial(i,4), 'k.')
        hold on
    end
    title(['Spatial = ', num2str(spatial_tested(t))])
    xlabel('Direction')
    ylabel('Spikes count during stimulus')
    
    clear spat64
   end
   
   
    
    suptitle({run_opt.data_set;['Cell ', num2str(cell_to_run)]});
    
    
    if ~exist(filepath)
        mkdir(filepath)
    end
    
    % set(gcf, 'renderer', 'opengl')
    
    export_fig([filepath, 'Cell_',num2str(cell_to_run)], '-pdf')
    % print(gcf,[filepath, 'Cell_',num2str(cell_to_run)], '-depsc', '-opengl');
    close(gcf);
end
