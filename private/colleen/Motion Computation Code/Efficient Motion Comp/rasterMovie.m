% make raster movie
clear toPlot
% DATA PARAMETERS
run_opt.load = true; % T/F

% run_opt.data_set = '2007-03-27-1';
run_opt.data_set = '2007-08-24-4';

run_opt.data_run = 6; % 12-19 for 2007-03-27, 2-11 for 2007-08-24, 13-17 for 2005-04-26
% CHANGE THIS
run_opt.config_num = 3; % 1-4 %Which type of stimulus to look at

run_opt.cell_type = 'Off parasol'; % on/off parasol, on/off midget

run_opt.velocity_exp = 190;

run_opt.cell_types = {'Off midget', 'Off parasol', 'On midget', 'On parasol'};

run_opt.auto_set = false; % T/F -- note: overwrites run_opt params

% NUMERICAL PARAMETERS
run_opt.tau = .01; % tuning parameter %0.1 next best
run_opt.tol = 1e-4;

% ANALYSES TO RUN
run_opt.downsample_spikes = false; % must run on bertha
run_opt.raster = true; % T/F
run_opt.rasterPerTrial = false; % T/F
run_opt.trial_estimate = false; % T/F

tic;

% Auto set parameters if flag set to true
if run_opt.auto_set
    [run_opt.cell_types, run_opt.velocity_lim, run_opt.config_num, run_opt.trial_estimate_start, run_opt.tol] =...
        auto_set_params(run_opt.data_set, run_opt.data_run);
end

% Load data
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
    
    % Get indices for specified cell type and order by RF position
    cell_indices1=get_cell_indices(datarun{1},{run_opt.cell_type});
    cell_indices2=get_cell_indices(datarun{2},{run_opt.cell_type});
    cell_x_pos = cellfun( @(X) X.mean(1), datarun{1}.vision.sta_fits); % x axis position of all STA cells
    [~, cell_sort_idx] = sort(cell_x_pos(cell_indices1)); % indicies of how to sort
    
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


% Plot one cell on all trials
if run_opt.raster %raster
    figure, set(gcf, 'Color','white')
    set(gca, 'nextplot','replacechildren', 'Visible','off');
    
    nFrames = length(cell_indices2);
    f = getframe(gca);
[f,map] = rgb2ind(f.cdata, 256, 'nodither');
mov = repmat(f, [1 1 1 nFrames]);
    for k = 1:length(cell_indices2)
        
        % Takes in start and stop time (0-0.7274)
        % Spikes of the cell with the lowest firing rate first
        % start time of each stimulus type trigger
        % Finds the spikes that happened on a cell from stimulus onset to end
        % Plot those spike times on the x axis versus the trial number on the y
        % axis
        % If tracking motion, the cell should respond to the bar at the same
        % time on every trial
        
        psth_r = psth_raster_noPlotting(start,stop,datarun{2}.spikes{cell_indices2(k)}',tr);
        Color = 'b.';
        y_scale =1;
        plot(psth_r(:,1),psth_r(:,2)*y_scale,Color);%, 'MarkerSize',10
        %     axis([mmin*1000 mmax*1000 0 length(tr)*y_scale]);
        
        % Title is the cell id according to vision and the mean firing rate
        title(sprintf('%d %.2f', datarun{2}.cell_ids(cell_indices2(k)), datarun{1}.vision.sta_fits{cell_indices1(k)}.mean(1) ))
        
        f = getframe(gca);
    mov(:,:,1,k) = rgb2ind(f.cdata, map, 'nodither');
    end
    close(gcf)
    
    %# create GIF and open
    imwrite(mov, map, 'myPeaks4.gif', 'DelayTime',0, 'LoopCount',inf)
    winopen('myPeaks4.gif')
    
end


% %# figure
% figure, set(gcf, 'Color','white')
% Z = peaks; surf(Z);  axis tight
% set(gca, 'nextplot','replacechildren', 'Visible','off');
% 
% %# preallocate
% nFrames = 20;
% f = getframe(gca);
% [f,map] = rgb2ind(f.cdata, 256, 'nodither');
% mov = repmat(f, [1 1 1 nFrames]);
% 
% %# create movie
% for k=1:nFrames
%     surf(sin(2*pi*k/20)*Z, Z)
%     f = getframe(gca);
%     mov(:,:,1,k) = rgb2ind(f.cdata, map, 'nodither');
% end
% close(gcf)
% 
% %# create GIF and open
% imwrite(mov, map, 'myPeaks4.gif', 'DelayTime',0, 'LoopCount',inf)
% winopen('myPeaks4.gif')