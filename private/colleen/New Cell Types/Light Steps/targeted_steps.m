clear

%%%%%%%%%%%%%%%%%%%%%%%%%% INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load the Data for the DS run
run_opt.data_set = '2015-11-09-7';
run_opt.data_run = 'd09-20/data012';
location = 'Data';

interval = 2; % number of seconds each stimulus was displayed for 

% Where to save the data
filepath= ['/Users/colleen/Desktop/Light Steps/', run_opt.data_set, '/', run_opt.data_run, '/'];

% You can give the cells as all of them (datarun.cell_ids) or give
% specific vision ids
% Find the cell to run by mapping a large cell EI from a white noise run
cells = [1636 451 4531 5748 5915]; % 'all' or give a vector of vision cell ids
map_order = [151 1636 451 4531 5748 5915];

%%%%%%%%%%%%%%%%%%%%%%%%%% END INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    if ~exist(filepath)
        mkdir(filepath)
    end

    
% datarun{1}.names.rrs_params_path=['/Volumes/Analysis/', run_opt.data_set, '/', run_opt.data_run, '/', run_opt.data_run(end-6:end), '.params'];
% datarun{1}.names.rrs_sta_path=['/Volumes/Analysis/', run_opt.data_set, '/', run_opt.data_run, '/', 'data022', '.sta'];

datarun.names.rrs_neurons_path=['/Volumes/Analysis/', run_opt.data_set, '/', run_opt.data_run, '/', run_opt.data_run(end-6:end), '.neurons'];
datarun.names.stimulus_path=['/Volumes/',location,'/', run_opt.data_set, '/Visual/s', run_opt.data_run(end-1:end)];


opt=struct('verbose',1,'load_sta', 1,'load_params',0,'load_neurons',1,'load_obvius_sta_fits',true);

datarun=load_data(datarun,opt);
% datarun=map_cell_types(datarun, struct('map',[1 2],'verbose',true));

% If the trigger interval in LabView is set long enought (~6 seconds for 5
% second stimuli), then this trigger_iti_thr should be fine.
datarun=load_stim(datarun,'correction_incomplet_run', 0, 'trigger_iti_thr', 0.1); % manually set threshold until got right number of triggers

if strcmp(cells, 'all')
    cells = datarun.cell_ids;
end

cell_index = get_cell_indices(datarun, cells);
for x = 1:length(cells)
    map_ind = find(map_order == cells(x)) - 1;
    trials_by_cell = [];
    for i = 1:length(datarun.stimulus.trials)
        if datarun.stimulus.trials(i).MAP == map_ind
            trials_by_cell = [trials_by_cell i];
        end
    end
    
    
    
    trial_ids = unique(datarun.stimulus.trial_list(trials_by_cell));
    step1_trials = find(datarun.stimulus.trial_list == trial_ids(1));
    step2_trials = find(datarun.stimulus.trial_list == trial_ids(2));
    
    step1_triggers = datarun.stimulus.triggers(step1_trials);
    step2_triggers = datarun.stimulus.triggers(step2_trials);
    
    spikes = datarun.spikes{cell_index(x)};
    cnt = 0;
    rasters1 = [];
        rasters2 = [];

    for trig = 1:length(step1_triggers)
        tmp=spikes(spikes>step1_triggers(trig) & spikes< step1_triggers(trig)+interval);
        rasters1=[rasters1, (tmp'-step1_triggers(trig))+interval*(trig-1)];
        cnt = cnt+1;
    end
%   L samples, Sampling rate = 1kHz. Spiketimes are hashed by the trial length.
    cnt = 0;
   for trig = 1:length(step2_triggers)
        tmp=spikes(spikes>step2_triggers(trig) & spikes< step2_triggers(trig)+interval);
        rasters2=[rasters2, (tmp'-step2_triggers(trig))+interval*(trig-1)];
        cnt = cnt+1;
    end
    
    
    
        figure;    
    if datarun.stimulus.trials(step1_trials(1)).RGB(1) > 0

        h= subplot(2,1,1);
        rasterplot(rasters1,25 ,interval, h)
        set(gca,'xlim', [0 ,interval])
        title('Step Up                                     Step Down')
        hold on
                y_lim = get(gca, 'ylim');

        plot(zeros(10,1), linspace(y_lim(1), y_lim(2) ,10)', 'g--')
        
        plot(ones(1,10), linspace(y_lim(1), y_lim(2) ,10), 'r--')
        
        g =subplot(2,1,2);
        rasterplot(rasters2,25 ,interval, g)
        set(gca,'xlim', [0 ,interval])
        title('Step Down                                    Step Up')
        hold on
                y_lim = get(gca, 'ylim');

        plot(zeros(1,10), linspace(y_lim(1), y_lim(2) ,10), 'r--')
        plot(ones(1,10), linspace(y_lim(1), y_lim(2) ,10), 'g--')
    else
        h = subplot(2,1,2);
        rasterplot(rasters1,25 ,interval,h)
        set(gca,'xlim', [0 ,interval]) 
        title('Step Down                                    Step Up')
        hold on
                y_lim = get(gca, 'ylim');

        plot(zeros(1,10), linspace(y_lim(1), y_lim(2) ,10), 'r--')
        plot(ones(1,10), linspace(y_lim(1), y_lim(2) ,10), 'g--')
        
        g = subplot(2,1,1)    ;   
        rasterplot(rasters2,25 ,interval, g)
        set(gca,'xlim', [0 ,interval])
        title('Step Up                                     Step Down')
        hold on
        y_lim = get(gca, 'ylim');
        plot(zeros(1,10), linspace(y_lim(1), y_lim(2) ,10), 'g--')
        plot(ones(1,10), linspace(y_lim(1), y_lim(2) ,10), 'r--')
        
    end
    suptitle({run_opt.data_set;[' Cell ' num2str(cells(x))]})
    print(gcf,'-dpdf',[filepath, 'Cell_',num2str(cells(x))]);
    
end




