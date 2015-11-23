clear

%%%%%%%%%%%%%%%%%%%%%%%%%% INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load the Data for the DS run
run_opt.data_set = '2015-10-06-5';
run_opt.data_run = 'data004-from-data002';
location = 'Data';

interval = 50; % number of seconds each stimulus was displayed for 

% Where to save the data
filepath= ['/Users/colleen/Desktop/Light Steps/Full Field', run_opt.data_set, '/', run_opt.data_run, '/'];

% You can give the cells as all of them (datarun{2}.cell_ids) or give
% specific vision ids
% Find the cell to run by mapping a large cell EI from a white noise run
cells = [3755]; % 'all' or give a vector of vision cell ids
% map_order = [1954 3107 3244 364 3679 4548 4850 4853 5057 5107 6136 6213 6406 6529 6646 6992 7043 723 7447];

%%%%%%%%%%%%%%%%%%%%%%%%%% END INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    if ~exist(filepath)
        mkdir(filepath)
    end

    
datarun{1}.names.rrs_params_path=['/Volumes/Analysis/', run_opt.data_set, '/', run_opt.data_run, '/', 'data004-from-data002', '.params'];
% datarun{1}.names.rrs_sta_path=['/Volumes/Analysis/', run_opt.data_set, '/', run_opt.data_run, '/', 'data022', '.sta'];

datarun{2}.names.rrs_neurons_path=['/Volumes/Analysis/', run_opt.data_set, '/', run_opt.data_run, '/', 'data004-from-data002', '.neurons'];

opt=struct('verbose',1,'load_sta', 1,'load_params',1,'load_neurons',1,'load_obvius_sta_fits',true);

datarun=load_data(datarun,opt);
datarun=map_cell_types(datarun, struct('map',[1 2],'verbose',true));

% If the trigger interval in LabView is set long enought (~6 seconds for 5
% second stimuli), then this trigger_iti_thr should be fine.
if strcmp(cells, 'all')
    cells = datarun{2}.cell_ids;
end

cell_index = get_cell_indices(datarun{1}, cells);
for x = 1:length(cells)
 
triggers = datarun{2}.triggers;
    on_triggers = triggers(1:2:end);
    off_triggers = triggers(2:2:end);
    
    spikes = datarun{2}.spikes{cell_index(x)};
    cnt = 0;
    rasters1 = [];
        rasters2 = [];

    for trig = 1:length(on_triggers)
        if trig == length(on_triggers)
           tmp=spikes(spikes>on_triggers(trig));

        else
        
        tmp=spikes(spikes>on_triggers(trig) & spikes< on_triggers(trig+1));
        end
        
        rasters1=[rasters1, (tmp'-on_triggers(trig))+interval*(trig-1)];
        cnt = cnt+1;
    end
%   L samples, Sampling rate = 1kHz. Spiketimes are hashed by the trial length.
%     cnt = 0;
%    for trig = 1:length(step2_triggers)
%         tmp=spikes(spikes>step2_triggers(trig) & spikes< step2_triggers(trig)+interval);
%         rasters2=[rasters2, (tmp'-step2_triggers(trig))+interval*(trig-1)];
%         cnt = cnt+1;
%     end
%     
    
    
        figure;    

        rasterplot(rasters1,10 ,interval)
        set(gca,'xlim', [0 ,interval])
        title('Step Up                                     Step Down')
        hold on
                y_lim = get(gca, 'ylim');

        plot(zeros(10,1), linspace(y_lim(1), y_lim(2) ,10)', 'g--')
        
        plot(10*ones(1,10), linspace(y_lim(1), y_lim(2) ,10), 'r--')
        
       
    suptitle({run_opt.data_set;[' Cell ' num2str(cells(x))]})
    print(gcf,'-dpdf',[filepath, 'Cell_',num2str(cells(x))]);
    
end




