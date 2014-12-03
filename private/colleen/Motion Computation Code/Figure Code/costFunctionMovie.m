% costFunctionMovie

% make raster movie
clear toPlot
% DATA PARAMETERS
run_opt.load = true; % T/F

run_opt.data_set = '2007-03-27-1';
% run_opt.data_set = '2007-08-24-4';

run_opt.data_run = 12; % 12-19 for 2007-03-27, 2-11 for 2007-08-24, 13-17 for 2005-04-26
% CHANGE THIS
run_opt.config_num = 14; % 1-4 %Which type of stimulus to look at

run_opt.cell_type = 'Off midget'; % on/off parasol, on/off midget

run_opt.velocity_exp = 96;

run_opt.cell_types = {'Off midget', 'Off parasol', 'On midget', 'On parasol'};

run_opt.auto_set = false; % T/F -- note: overwrites run_opt params

direction = 'right';
% NUMERICAL PARAMETERS
run_opt.tau = .01; % tuning parameter %0.1 next best
run_opt.tol = 1e-4;

% ANALYSES TO RUN
run_opt.downsample_spikes = false; % must run on bertha
run_opt.raster = false; % T/F
run_opt.rasterPerTrial = false; % T/F
run_opt.trial_estimate = true; % T/F

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
    
    writerObj = VideoWriter('Results/resultsColleen/raster_alltrials.avi');
    writerObj.FrameRate = 8;
    
    open(writerObj);
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
        xlim([200 800]);
        
        %     axis([mmin*1000 mmax*1000 0 length(tr)*y_scale]);
        
        % Title is the cell id according to vision and the mean firing rate
        title(sprintf('%d %.2f', datarun{2}.cell_ids(cell_indices2(k)), datarun{1}.vision.sta_fits{cell_indices1(k)}.mean(1) ))
        
        frame = getframe;
        writeVideo(writerObj,frame);
        
        
    end
    close(gcf)
    close(writerObj);
    
    
end

if run_opt.rasterPerTrial
    toPlot = cell(1,length(t));
    figure, set(gcf, 'Color','white')
    set(gca, 'nextplot','replacechildren', 'Visible','off');
    
    writerObj = VideoWriter('Results/resultsColleen/raster_allcells.avi');
    writerObj.FrameRate = 8;
    
    open(writerObj);
    % Takes in start and stop time (0-0.7274)
    % Spikes of the cell with the lowest firing rate first
    % start time of each stimulus type 2 trigger
    % Finds the spikes that happened on a cell from stimulus onset to end
    % Plot those spike times on the x axis versus the trial number on the y
    % axis
    % If tracking motion, the cell should respond to the bar at the same
    % time on every trial
    for counter = 1:length(t)
        psth_r = psth_raster_noPlotting(start,stop,datarun{2}.spikes{cell_indices2(counter)}',tr);
        posThisCell = datarun{1}.vision.sta_fits{cell_indices1(counter)}.mean(1);
        
        posFarthestCell = datarun{1}.vision.sta_fits{cell_indices1(1)}.mean(1);
        
        
        cellNumber = datarun{2}.cell_ids(cell_indices2(counter));
        % Title is the cell id according to vision and the mean firing rate
        %          [psth, bins] = get_psth(datarun{2}.spikes{cell_indices2(counter)}, tr, 'plot_hist', true)
        for trialNum = 1:length(t)
            [x,y] = find(psth_r == trialNum-1);
            toPlot{trialNum}= [toPlot{trialNum}; [psth_r(x,1), repmat(cellNumber, length(x),1), repmat(posThisCell, length(x),1)]];
        end
    end
    
    
    for k = 1:length(t)
        y_scale = 1;
        Color = ['k', '.'];
        plot(toPlot{k}(:,1),toPlot{k}(:,3)*y_scale,Color);
        xlim([0 800])
        title({run_opt.cell_type, [run_opt.data_set, ' Run ', num2str(run_opt.data_run)],'Bright Bars Moving Right', sprintf(' Trial Number %d',  k)})
        xlabel('time (ms)');
        ylabel('Cell''s centroid distance from reference');
        frame = getframe;
        writeVideo(writerObj,frame);
    end
    close(gcf)
    close(writerObj);
    
end

if run_opt.trial_estimate
    options = optimset('Display', 'iter', 'TolFun', run_opt.tol , 'MaxFunEvals', 60, 'LargeScale', 'off');
    toPlot = cell(1,length(t));
    figure, set(gcf, 'Color','white')
    set(gca, 'nextplot','replacechildren', 'Visible','off');
    
    writerObj = VideoWriter('Results/resultsColleen/costFunction_OFFM.avi');
    writerObj.FrameRate = 8;
    
    open(writerObj);
    
    spikes = datarun{2}.spikes;
    
    %Prior is +/-25% of expected value
    velocity = linspace(0.75*run_opt.velocity_exp, 1.25*run_opt.velocity_exp, 50);
    %     velocity = 100:5:400;
    strsig1 = zeros(1,length(velocity));
    
    % Run coarse error function to initialize velocity
    for i =1:length(tr)
        parfor j = 1:length(velocity)
            v = velocity(j);
            [strsig1(j)] = pop_motion_signal_colleen(v, spikes, cell_indices1, cell_indices2, cell_x_pos, tr(i), stop, run_opt.tau, run_opt.tol, datarun, direction);
        end
        plot(velocity, strsig1)
        axis off
        i
        [x1,y1] = min(strsig1);
        
        % Initialize minimization
        frame = getframe;
        frame.cdata = frame.cdata(1:343, 1:435, :);
        writeVideo(writerObj,frame);
    end
    
    close(writerObj);
    
end




