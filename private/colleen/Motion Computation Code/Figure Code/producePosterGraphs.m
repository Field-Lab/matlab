clear
global savedVariables
close all
run_opt.load = true; % T/F

% run_opt.data_set = '2007-03-27-1';
run_opt.data_set = '2007-08-24-4';

% run_opt.data_set = '2007-08-24-4';
run_opt.data_run = 7; % 12-19 for 2007-03-27, 2-11 for 2007-08-24, 13-17 for 2005-04-26
run_opt.config_num = 3; % 1-4 %Which type of stimulus to look at


% stim categories not consistant
%1: dark bar, x_delta= 8
%2 dark bar, x_delta = -8
%3 light bar, x_delta = 8
%4 light bar, x_delta = -8

% Change this to change type of cell you are interested in

run_opt.cell_type = 'On parasol'; % on/off parasol, on/off midget

run_opt.cell_types = {'Off midget', 'Off parasol', 'On midget', 'On parasol'};
run_opt.auto_set = false; % T/F -- note: overwrites run_opt params

% NUMERICAL PARAMETERS
run_opt.tau = .01; % tuning parameter
run_opt.tol = 1e-3;

run_opt.trial_estimate_start = 150;

run_opt.velocity_lim = 150; % >0

% ANALYSES TO RUN
run_opt.downsample_spikes = false; % must run on bertha
run_opt.raster = true; % T/F
run_opt.rasterPerTrial = true; % T/F
run_opt.trial_estimate = false; % T/F

speed =0.09;

tic;

% Auto set parameters
if run_opt.auto_set
    [run_opt.cell_types, run_opt.velocity_lim, run_opt.config_num, run_opt.trial_estimate_start, run_opt.tol] =...
        auto_set_params(run_opt.data_set, run_opt.data_run);
end

% Load data fresh
if run_opt.load
    
    clear datarun tr
    % datarun{1} has vision info (sta fits)
    % datarun{2} has cell_ids, spikes, triggers
    if strcmp(run_opt.data_set, '2007-03-27-1')
        datarun{1}.names.rrs_params_path='/Volumes/Analysis/2007-03-27-1/data011-nwpca/data011-nwpca.params';
        datarun{2}.names.rrs_neurons_path=sprintf('/Volumes/Analysis/2007-03-27-1/data%03d-from-data011-nwpca/data%03d-from-data011-nwpca.neurons', run_opt.data_run, run_opt.data_run);
        datarun{2}.names.stimulus_path=sprintf('/Volumes/Analysis/2007-03-27-1/stimuli/s%02d', run_opt.data_run);
    elseif strcmp(run_opt.data_set, '2007-08-24-4')
        datarun{1}.names.rrs_params_path='/Volumes/Analysis/2007-08-24-4/data001-nwpca/data001-nwpca.params';
        datarun{2}.names.rrs_neurons_path=sprintf('/Volumes/Analysis/2007-08-24-4/data%03d-from-data001-nwpca/data%03d-from-data001-nwpca.neurons', run_opt.data_run, run_opt.data_run);
        datarun{2}.names.stimulus_path=sprintf('/Volumes/Analysis/2007-08-24-4/Stimuli/s%02d', run_opt.data_run);
    elseif strcmp(run_opt.data_set, '2005-04-26-0')
        datarun{1}.names.rrs_params_path='/Volumes/Analysis/2005-04-26-0/data';
    end
    opt=struct('verbose',1,'load_params',1,'load_neurons',1,'load_obvius_sta_fits',true);
    datarun=load_data(datarun,opt);
    datarun=map_cell_types(datarun, struct('map',[1 2],'verbose',true));
    datarun{2}=load_stim(datarun{2},'correction_incomplet_run', 0);
    
end



% Gets the indicies used by vision of the particular cell type
if run_opt.raster || run_opt.trial_estimate || run_opt.rasterPerTrial
    
    
    % ex. in vision, on midgets have IDs 17, 61, 110 etc
    % cell_indicies1 is the index of the cell_ids array containing vision
    % ids that corresponds to the on midget cells
    % Get indices for specified cell type and order by RF position
    cell_indices1=get_cell_indices(datarun{1},{run_opt.cell_type});
    cell_indices2=get_cell_indices(datarun{2},{run_opt.cell_type});
    true
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
    
    
end


if run_opt.rasterPerTrial %raster
    switch run_opt.config_num % 1-4 %Which type of stimulus to look at
        case 1
            strTit = 'Bright Bar Moving Right';
        case 2
            strTit = 'Bright Bar Moving Right';
        case 3
            strTit = 'Bright Bar Moving Right';
        case 4
            strTit = 'Bright Bar Moving Right';
    end
    
    toPlot = cell(1,length(t));
    
    % Takes in start and stop time (0-0.7274)
    % Spikes of the cell with the lowest firing rate first
    % start time of each stimulus type 2 trigger
    % Finds the spikes that happened on a cell from stimulus onset to end
    % Plot those spike times on the x axis versus the trial number on the y
    % axis
    % If tracking motion, the cell should respond to the bar at the same
    % time on every trial
    for counter = 2:length(cell_indices2)
        psth_r = psth_raster_noPlotting(start,stop,datarun{2}.spikes{cell_indices2(counter)}',tr);
        posThisCell = datarun{1}.vision.sta_fits{cell_indices1(counter)}.mean(1);
        
        posFarthestCell = datarun{1}.vision.sta_fits{cell_indices1(1)}.mean(1);
        
        
        cellNumber = datarun{2}.cell_ids(cell_indices2(counter));
        % Title is the cell id according to vision and the mean firing rate
        %          [psth, bins] = get_psth(datarun{2}.spikes{cell_indices2(counter)}, tr, 'plot_hist', true)
        for trialNum = 1:length(t)
            [x,y] = find(psth_r == trialNum-1);
            if run_opt.config_num == 1 || run_opt.config_num == 3
                %                         toPlot{trialNum}= [toPlot{trialNum}; [psth_r(x,1)-repmat(posThisCell - posFarthestCell, length(x),1)/speed, repmat(cellNumber, length(x),1), repmat(posThisCell, length(x),1)]];
                toPlot{trialNum}= [toPlot{trialNum}; [psth_r(x,1), repmat(cellNumber, length(x),1), repmat(posThisCell, length(x),1)]];
                
            else
                %                        toPlot{trialNum}= [toPlot{trialNum}; [psth_r(x,1)-repmat(-posThisCell + posFarthestCell, length(x),1)/speed, repmat(cellNumber, length(x),1), repmat(posThisCell, length(x),1)]];
                toPlot{trialNum}= [toPlot{trialNum}; [psth_r(x,1), repmat(cellNumber, length(x),1), repmat(posThisCell, length(x),1)]];
                
            end
            
        end
    end
    
    k=1; kmin=1; kmax=length(t); hk=1;
    
    hFig = figure(1);
    set(hFig, 'Position', [680   678   300   420])
    k=1;
    y_scale = 1;
    Color = ['k', '.'];
    plot(toPlot{k}(:,1),toPlot{k}(:,3)*y_scale,Color, 'MarkerSize', 15);%, 'MarkerSize',10
    xlim([0 775])
    title({'ON Parasol'; 'Faster Speed';strTit })
    xlabel('time (ms)');
    ylabel('Cell''s centroid distance from reference');
    ylim([15 50])
    set(gcf,'color','w');
end
%%%%%%%%%%%% graph 2



run_opt.load = true; % T/F

% run_opt.data_set = '2007-03-27-1';
run_opt.data_set = '2007-08-24-4';

% run_opt.data_set = '2007-08-24-4';
run_opt.data_run = 10; % 12-19 for 2007-03-27, 2-11 for 2007-08-24, 13-17 for 2005-04-26
run_opt.config_num = 2; % 1-4 %Which type of stimulus to look at


% stim categories not consistant
%1: dark bar, x_delta= 8
%2 dark bar, x_delta = -8
%3 light bar, x_delta = 8
%4 light bar, x_delta = -8

% Change this to change type of cell you are interested in

run_opt.cell_type = 'On parasol'; % on/off parasol, on/off midget

run_opt.cell_types = {'Off midget', 'Off parasol', 'On midget', 'On parasol'};
run_opt.auto_set = false; % T/F -- note: overwrites run_opt params

% NUMERICAL PARAMETERS
run_opt.tau = .01; % tuning parameter
run_opt.tol = 1e-3;

run_opt.trial_estimate_start = 150;

run_opt.velocity_lim = 150; % >0

% ANALYSES TO RUN
run_opt.downsample_spikes = false; % must run on bertha
run_opt.raster = true; % T/F
run_opt.rasterPerTrial = true; % T/F
run_opt.trial_estimate = false; % T/F

speed =0.09;

tic;

% Auto set parameters
if run_opt.auto_set
    [run_opt.cell_types, run_opt.velocity_lim, run_opt.config_num, run_opt.trial_estimate_start, run_opt.tol] =...
        auto_set_params(run_opt.data_set, run_opt.data_run);
end

% Load data fresh
if run_opt.load
    
    clear datarun tr
    % datarun{1} has vision info (sta fits)
    % datarun{2} has cell_ids, spikes, triggers
    if strcmp(run_opt.data_set, '2007-03-27-1')
        datarun{1}.names.rrs_params_path='/Volumes/Analysis/2007-03-27-1/data011-nwpca/data011-nwpca.params';
        datarun{2}.names.rrs_neurons_path=sprintf('/Volumes/Analysis/2007-03-27-1/data%03d-from-data011-nwpca/data%03d-from-data011-nwpca.neurons', run_opt.data_run, run_opt.data_run);
        datarun{2}.names.stimulus_path=sprintf('/Volumes/Analysis/2007-03-27-1/stimuli/s%02d', run_opt.data_run);
    elseif strcmp(run_opt.data_set, '2007-08-24-4')
        datarun{1}.names.rrs_params_path='/Volumes/Analysis/2007-08-24-4/data001-nwpca/data001-nwpca.params';
        datarun{2}.names.rrs_neurons_path=sprintf('/Volumes/Analysis/2007-08-24-4/data%03d-from-data001-nwpca/data%03d-from-data001-nwpca.neurons', run_opt.data_run, run_opt.data_run);
        datarun{2}.names.stimulus_path=sprintf('/Volumes/Analysis/2007-08-24-4/Stimuli/s%02d', run_opt.data_run);
    elseif strcmp(run_opt.data_set, '2005-04-26-0')
        datarun{1}.names.rrs_params_path='/Volumes/Analysis/2005-04-26-0/data';
    end
    opt=struct('verbose',1,'load_params',1,'load_neurons',1,'load_obvius_sta_fits',true);
    datarun=load_data(datarun,opt);
    datarun=map_cell_types(datarun, struct('map',[1 2],'verbose',true));
    datarun{2}=load_stim(datarun{2},'correction_incomplet_run', 0);
    
end



% Gets the indicies used by vision of the particular cell type
if run_opt.raster || run_opt.trial_estimate || run_opt.rasterPerTrial
    
    
    % ex. in vision, on midgets have IDs 17, 61, 110 etc
    % cell_indicies1 is the index of the cell_ids array containing vision
    % ids that corresponds to the on midget cells
    % Get indices for specified cell type and order by RF position
    cell_indices1=get_cell_indices(datarun{1},{run_opt.cell_type});
    cell_indices2=get_cell_indices(datarun{2},{run_opt.cell_type});
    true
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
    
    
end


if run_opt.rasterPerTrial %raster
    switch run_opt.config_num % 1-4 %Which type of stimulus to look at
        case 1
            strTit = 'Bright Bar Moving Right';
        case 2
            strTit = 'Bright Bar Moving Right';
        case 3
            strTit = 'Bright Bar Moving Right';
        case 4
            strTit = 'Bright Bar Moving Right';
    end
    
    toPlot = cell(1,length(t));
    
    % Takes in start and stop time (0-0.7274)
    % Spikes of the cell with the lowest firing rate first
    % start time of each stimulus type 2 trigger
    % Finds the spikes that happened on a cell from stimulus onset to end
    % Plot those spike times on the x axis versus the trial number on the y
    % axis
    % If tracking motion, the cell should respond to the bar at the same
    % time on every trial
    for counter = 2:length(cell_indices2)
        psth_r = psth_raster_noPlotting(start,stop,datarun{2}.spikes{cell_indices2(counter)}',tr);
        posThisCell = datarun{1}.vision.sta_fits{cell_indices1(counter)}.mean(1);
        
        posFarthestCell = datarun{1}.vision.sta_fits{cell_indices1(1)}.mean(1);
        
        
        cellNumber = datarun{2}.cell_ids(cell_indices2(counter));
        % Title is the cell id according to vision and the mean firing rate
        %          [psth, bins] = get_psth(datarun{2}.spikes{cell_indices2(counter)}, tr, 'plot_hist', true)
        for trialNum = 1:length(t)
            [x,y] = find(psth_r == trialNum-1);
            if run_opt.config_num == 1 || run_opt.config_num == 3
                %                         toPlot{trialNum}= [toPlot{trialNum}; [psth_r(x,1)-repmat(posThisCell - posFarthestCell, length(x),1)/speed, repmat(cellNumber, length(x),1), repmat(posThisCell, length(x),1)]];
                toPlot{trialNum}= [toPlot{trialNum}; [psth_r(x,1), repmat(cellNumber, length(x),1), repmat(posThisCell, length(x),1)]];
                
            else
                %                        toPlot{trialNum}= [toPlot{trialNum}; [psth_r(x,1)-repmat(-posThisCell + posFarthestCell, length(x),1)/speed, repmat(cellNumber, length(x),1), repmat(posThisCell, length(x),1)]];
                toPlot{trialNum}= [toPlot{trialNum}; [psth_r(x,1), repmat(cellNumber, length(x),1), repmat(posThisCell, length(x),1)]];
                
            end
            
        end
    end
    
    k=1; kmin=1; kmax=length(t); hk=1;
    
    
    k=1;
    y_scale = 1;
    Color = ['k', '.'];
    hFig = figure(2);
    set(hFig, 'Position', [680   678   560   420])
    plot(toPlot{k}(:,1),toPlot{k}(:,3)*y_scale,Color, 'MarkerSize', 15);%, 'MarkerSize',10
    xlim([0 1500])
    title({'ON Parasol'; 'Slower Speed'; 'Bright Bar Moving Right'})
    xlabel('time (ms)');
    ylabel('Cell''s centroid distance from reference');
    ylim([15 50])
    set(gcf,'color','w');
    
    
end



%% Error function


run_opt.load = true; % T/F

% run_opt.data_set = '2007-03-27-1';
run_opt.data_set = '2007-08-24-4';

run_opt.data_run = 4; % 12-19 for 2007-03-27, 2-11 for 2007-08-24, 13-17 for 2005-04-26
run_opt.config_num = 1; % 1-4 %Which type of stimulus to look at

run_opt.cell_type = 'On midget'; % on/off parasol, on/off midget
run_opt.direction = 'right';
run_opt.cell_types = {'Off midget', 'Off parasol', 'On midget', 'On parasol'};
run_opt.auto_set = false; % T/F -- note: overwrites run_opt params

% NUMERICAL PARAMETERS
run_opt.tau = .01; % tuning parameter %0.1 next best
run_opt.tol = 1e-3;
% ANALYSES TO RUN
run_opt.downsample_spikes = false; % must run on bertha
run_opt.raster = false; % T/F
run_opt.rasterPerTrial = false; % T/F
run_opt.trial_estimate = true; % T/F

  if strcmp(run_opt.data_set, '2007-03-27-1')
        run_opt.stx = 8;
    elseif strcmp(run_opt.data_set, '2007-08-24-4')
        run_opt.stx = 10;
    else
        run_opt.stx = 10;
  end
    
% Load data fresh
if run_opt.load
    
    clear datarun tr
    % datarun{1} has vision info (sta fits)
    % datarun{2} has cell_ids, spikes, triggers
    if strcmp(run_opt.data_set, '2007-03-27-1')
        datarun{1}.names.rrs_params_path='/Volumes/Analysis/2007-03-27-1/data011-nwpca/data011-nwpca.params';
        datarun{2}.names.rrs_neurons_path=sprintf('/Volumes/Analysis/2007-03-27-1/data%03d-from-data011-nwpca/data%03d-from-data011-nwpca.neurons', run_opt.data_run, run_opt.data_run);
        datarun{2}.names.stimulus_path=sprintf('/Volumes/Analysis/2007-03-27-1/stimuli/s%02d', run_opt.data_run);
    elseif strcmp(run_opt.data_set, '2007-08-24-4')
        datarun{1}.names.rrs_params_path='/Volumes/Analysis/2007-08-24-4/data001-nwpca/data001-nwpca.params';
        datarun{2}.names.rrs_neurons_path=sprintf('/Volumes/Analysis/2007-08-24-4/data%03d-from-data001-nwpca/data%03d-from-data001-nwpca.neurons', run_opt.data_run, run_opt.data_run);
        datarun{2}.names.stimulus_path=sprintf('/Volumes/Analysis/2007-08-24-4/Stimuli/s%02d', run_opt.data_run);
    elseif strcmp(run_opt.data_set, '2005-04-26-0')
        datarun{1}.names.rrs_params_path='/Volumes/Analysis/2005-04-26-0/data';
    end
    opt=struct('verbose',1,'load_params',1,'load_neurons',1,'load_obvius_sta_fits',true);
    datarun=load_data(datarun,opt);
    datarun=map_cell_types(datarun, struct('map',[1 2],'verbose',true));
    datarun{2}=load_stim(datarun{2},'correction_incomplet_run', 0);
    
end


% ex. in vision, on midgets have IDs 17, 61, 110 etc
% cell_indicies1 is the index of the cell_ids array containing vision
% ids that corresponds to the on midget cells
% Get indices for specified cell type and order by RF position
cell_indices1=get_cell_indices(datarun{1},{run_opt.cell_type});
cell_indices2=get_cell_indices(datarun{2},{run_opt.cell_type});
true
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


spikes = datarun{2}.spikes;


velocity= [1:2:400];
strsig1 = zeros(1,length(velocity));

parfor j = 1:length(velocity)
    v = velocity(j);
    [strsig1(j)] = -pop_motion_signal(v, spikes, cell_indices1, cell_indices2, cell_x_pos, tr(23), stop, run_opt.tau, run_opt.tol, datarun, run_opt.direction, run_opt.stx);
end

[x1,y1] = min(strsig1);

run_opt.trial_estimate_start = velocity(y1);


hFig = figure;
set(hFig, 'Position', [680 680 600 435])
plot(velocity*5/225*10, -strsig1, 'k', 'LineWidth', 2')
s = sprintf('%c', char(176));
xlabel(['Speed Estimates (', s, '/sec)'])
ylabel('N: Motion Signal (arb. units)')
set(gca, 'ytick', []);
set(gca, 'yticklabel', []);
title({'ON Midget'; 'Motion Signal';'Bright Bar Moving Right'});
hold on
plot(velocity(y1)*5/225*10, -strsig1(y1), 'ro', 'MarkerFaceColor', 'r')
hold on
y_bottom = get(gca, 'ylim');
xlim([0 45])
plot(96*5/225*10,y_bottom(1), 'bo', 'MarkerFaceColor', 'b')

hold on
plot(19.65,41.37, 'go', 'MarkerFaceColor', 'g')
legend('Error Function', 'Estimated Speed', 'True Speed', 'Estimated Speed with Prior')

plot(velocity(y1)*5/225*10*ones(100,1), linspace(y_bottom(1), -strsig1(y1), 100),'--r', 'LineWidth', 2);
plot(19.65*ones(100,1), linspace(y_bottom(1), 41.37, 100),'--g', 'LineWidth', 2);

set(gcf,'color','w');


%%

load('/Users/vision/Desktop/GitHub code repository/private/colleen/Results/resultsColleen/2007-08-24-4/BrightRight/On parasol_data_run_10_config_2.mat')

figure;
histfit(estimates*5/225*10, 10)
% xlim([10.4 11.07]);
s = sprintf('%c', char(176));
xlabel(['Speed Estimates (', s, '/sec)'])
ylabel('Counts')
title({'ON Parasol'; 'Histogram of Speed Estimates'; 'Bright Bar Moving Right'});
children = get(gca,'children');
set(children(2),'FaceColor',[1 1 1],'EdgeColor','k');
set(gcf,'color','w');


load('/Users/vision/Desktop/GitHub code repository/private/colleen/Results/resultsColleen/2007-08-24-4/BrightRight/Off parasol_data_run_06_config_3.mat')
figure;
histfit(estimates*5/225*10, 10)
% xlim([10.4 11.07]);
s = sprintf('%c', char(176));
xlabel(['Speed Estimates (', s, '/sec)'])
ylabel('Counts')
title({'OFF Parasol'; 'Histogram of Speed Estimates'; 'Bright Bar Moving Right'});
children = get(gca,'children');
set(children(2),'FaceColor',[1 1 1],'EdgeColor','k');
set(gcf,'color','w');

load('/Users/vision/Desktop/GitHub code repository/private/colleen/colleen/Results/resultsColleen/2007-08-24-4/BrightRight/On midget_data_run_02_config_9.mat')
figure;
histfit(estimates*5/225*10, 10)
% xlim([10.4 11.07]);
s = sprintf('%c', char(176));
xlabel(['Speed Estimates (', s, '/sec)'])
ylabel('Counts')
title({'ON Midget'; 'Histogram of Speed Estimates'; 'Bright Bar Moving Right'});
children = get(gca,'children');
set(children(2),'FaceColor',[1 1 1],'EdgeColor','k');
set(gcf,'color','w');

