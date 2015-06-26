% ds_cells
clear


% Load the Data for the DS run
run_opt.data_set = '2015-05-27-5';
run_opt.data_run = 3; % 12-19 for 2007-03-27, 2-11 for 2007-08-24, 13-17 for 2005-04-26
if strcmp(run_opt.data_set, '2015-05-27-5')
    datarun{1}.names.rrs_params_path='/Volumes/Analysis/2015-05-27-5/data003/data003.params';
    datarun{2}.names.rrs_neurons_path=('/Volumes/Analysis/2015-05-27-5/data003/data003.neurons');
    datarun{2}.names.stimulus_path=('/Volumes/Data/2015-05-27-5/Visual/s03');
end
opt=struct('verbose',1,'load_params',1,'load_neurons',1,'load_obvius_sta_fits',true);
datarun=load_data(datarun,opt);
datarun=map_cell_types(datarun, struct('map',[1 2],'verbose',true));
datarun{2}=load_stim(datarun{2},'correction_incomplet_run', 0, 'trigger_iti_thr', 0.1); % manually set threshold until got right number of triggers
% Something went wrong with the first two stimuli
% datarun{2}.stimulus.triggers= [datarun{2}.stimulus.triggers(1:33) 187.0086, datarun{2}.stimulus.triggers(34:70) 396.0097 datarun{2}.stimulus.triggers(71:72) 412.5063 datarun{2}.stimulus.triggers(73:72); % Remove the first trigger
% datarun{2}.stimulus.trials = datarun{2}.stimulus.trials(3:end); % Remove the first two trials

interval = 5; % each stimulus was displayed for 5 seconds
spatial_tested = fliplr(datarun{2}.stimulus.params.SPATIAL_PERIOD);
temporal_tested = fliplr(datarun{2}.stimulus.params.TEMPORAL_PERIOD);
% Where to save the data
filepath= ['/Users/colleen/Desktop/DS/', run_opt.data_set, '/data003/'];


% Find the cell to run by mapping a large cell EI from a white noise run
cells = datarun{2}.cell_ids;
% cells = 513;

% get spikes per bin
    
    % These triggers depend on the Labview trigger. 7.5 seconds prevents
    % triggers from being missed on the large spatial scales
    tr=datarun{2}.stimulus.triggers; % all start triggers
    run_opt.config_num = 1;
    
for x = 1:length(cells)
    x
    cell_to_run = cells(x);
    cell_indices1=get_cell_indices(datarun{1}, [cell_to_run]);
    cell_indices2=get_cell_indices(datarun{2},[cell_to_run]);
    
    
    
    
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
