% DATA PARAMETERS
run_opt.load = true; % T/F
% run_opt.data_set = '2007-03-27-1';
run_opt.data_set = '2007-08-24-4';
run_opt.data_run = 7; % 12-19 for 2007-03-27, 2-11 for 2007-08-24, 13-17 for 2005-04-26
run_opt.config_num = 2; % 1-4 %Which type of stimulus to look at 
%1: dark bar, x_delta = 8
%2 dark bar, x_delta = -8
%3 light bar, x_delta = 8
%4 light bar, x_delta = -8

% Change this to change type of cell you are interested in
run_opt.cell_type = 'On midget'; % on/off parasol, on/off midget
run_opt.cell_types = {'Off midget', 'Off parasol', 'On midget', 'On parasol'};
run_opt.auto_set = false; % T/F -- note: overwrites run_opt params

% NUMERICAL PARAMETERS
run_opt.tau = .01; % tuning parameter
run_opt.tol = 1e-3;
run_opt.trial_estimate_start = 120;
run_opt.velocity_lim = 150; % >0

% ANALYSES TO RUN
run_opt.downsample_spikes = false; % must run on bertha
run_opt.raster = false; % T/F
run_opt.rasterPerTrial = true; % T/F
run_opt.trial_estimate = false; % T/F


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
    
    cell_x_pos = cellfun( @(X) X.mean(1), datarun{1}.vision.sta_fits); % mean firing rate of all STA cells
    [~, cell_sort_idx] = sort(cell_x_pos(cell_indices1)); % mean firing rate of only on midget cells, indexes of how to sort
    
    %cell_indices sorted by their 2nd entry of sta mean 
    cell_indices1 = cell_indices1(cell_sort_idx); % cell_indices1 is now indexes in order from lowest to highest firing rate
    cell_indices2 = cell_indices2(cell_sort_idx);
    
    % Find trial start and stop times
    start = 0;
    stop = mean(datarun{2}.triggers(2:2:end) - datarun{2}.triggers(1:2:end));
    tr=datarun{2}.triggers(1:2:end); % all start triggers
    t=find(datarun{2}.stimulus.trial_list==run_opt.config_num); %find the times when all the stimulus type 2 starts
    tr=tr(t);
    
end

% downsample spikes
if run_opt.downsample_spikes
    
    onp_indices = get_cell_indices(datarun{2}, {'On parasol'});
    numspikes1 = zeros(length(onp_indices),1);
    for i=1:length(onp_indices)
        numspikes1(i) = length(datarun{2}.spikes{onp_indices(i)});
    end
    
    onm_indices = get_cell_indices(datarun{2}, {'On midget'});
    numspikes2 = zeros(length(onm_indices),1);
    for i=1:length(onm_indices)
        numspikes2(i) = length(datarun{2}.spikes{onm_indices(i)});
    end
    
    n1 = mean(numspikes1);
    n2 = mean(numspikes2);
    
    if strcmp(run_opt.data_set, '2007-03-27-1')
        if n1>n2
            if strcmp(run_opt.cell_type, 'On parasol')
                for i=cell_indices2
                    % fractional downsampling:
                    datarun{2}.spikes{i}=datasample(datarun{2}.spikes{i},round(length(datarun{2}.spikes{i})*n2/n1),'Replace',false);
                end
            else
                run_opt.trial_estimate = false; % if no downsampling, don't calculate estimates
            end
        elseif n2>n1
            if strcmp(run_opt.cell_type, 'On midget')
                for i=cell_indices2
                    % fractional downsampling
                    datarun{2}.spikes{i}=datasample(datarun{2}.spikes{i},round(length(datarun{2}.spikes{i})*n1/n2),'Replace',false);                    
                end
            else
                run_opt.trial_estimate = false; % if no downsampling, don't calculate estimates
            end
        end
    elseif strcmp(run_opt.data_set, '2007-08-24-4')
        if n1>n2
            if strcmp(run_opt.cell_type, 'On parasol')
                for i=cell_indices2
                    % fractional downsampling
                    datarun{2}.spikes{i}=datasample(datarun{2}.spikes{i},round(length(datarun{2}.spikes{i})*n2/n1),'Replace',false);
                end
            else
                run_opt.trial_estimate = false; % if no downsampling, don't calculate estimates
            end
        elseif n2>n1    
            if strcmp(run_opt.cell_type, 'On midget')
                for i=cell_indices2
                    % fractional downsampling
                    datarun{2}.spikes{i}=datasample(datarun{2}.spikes{i},round(length(datarun{2}.spikes{i})*n1/n2),'Replace',false);
                end
            else
                run_opt.trial_estimate = false; % if no downsampling, don't calculate estimates
            end
        end
    end
end

if run_opt.raster %raster

    k=1; kmin=1; kmax=length(cell_indices2); hk=loop_slider(k,kmin,kmax);

    while k
        if ~ishandle(hk)
            break % script breaks until figure is closed
        end
        k=round(get(hk,'Value')); 
    % Takes in start and stop time (0-0.7274)
    % Spikes of the cell with the lowest firing rate first
    % start time of each stimulus type 2 trigger
    % Finds the spikes that happened on a cell from stimulus onset to end
    % Plot those spike times on the x axis versus the trial number on the y
    % axis
    % If tracking motion, the cell should respond to the bar at the same
    % time on every trial
        psth_r = psth_raster(start,stop,datarun{2}.spikes{cell_indices2(k)}',tr);
        
        % Title is the cell id according to vision and the mean firing rate
        title(sprintf('%d %.2f', datarun{2}.cell_ids(cell_indices2(k)), datarun{1}.vision.sta_fits{cell_indices1(k)}.mean(1) ))

        uiwait;
    end
end

if run_opt.rasterPerTrial %raster
    switch run_opt.config_num % 1-4 %Which type of stimulus to look at
        case 1
            strTit = 'Dark Bar Moving Left';
        case 2
            strTit = 'Dark Bar Moving Right';
        case 3
            strTit = 'Bright Bar Moving Left';
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
    for counter = 1:length(cell_indices2)
        psth_r = psth_raster_noPlotting(start,stop,datarun{2}.spikes{cell_indices2(counter)}',tr);
        
        % Title is the cell id according to vision and the mean firing rate
        for trialNum = 1:length(t)
            [x,y] = find(psth_r == trialNum-1);
            toPlot{trialNum}= [toPlot{trialNum}; [psth_r(x,1), repmat(datarun{2}.cell_ids(cell_indices2(counter)), length(x),1), repmat(datarun{1}.vision.sta_fits{cell_indices1(counter)}.mean(1), length(x),1)]];
        end
    end
    
     k=1; kmin=1; kmax=length(t); hk=loop_slider(k,kmin,kmax);

   while k
        if ~ishandle(hk)
            break % script breaks until figure is closed
        end
        k=round(get(hk,'Value')); 
    plot(toPlot{k}(:,1),toPlot{k}(:,3)*y_scale,Color);%, 'MarkerSize',10
    
    title({run_opt.cell_type, [run_opt.data_set, ' Run ', num2str(run_opt.data_run)],strTit, sprintf(' Trial Number %d',  k)})
 xlabel('time (ms)');
 ylabel('Centroid distance from reference');
    uiwait;
    end
    
end



if run_opt.trial_estimate
    
    % start parallel pool
    poolobj = parpool;
    
    options = optimset('Display', 'iter', 'TolFun', run_opt.tol , 'MaxFunEvals', 30, 'LargeScale', 'off');
    estimates = zeros(size(tr));
    spikes = datarun{2}.spikes;
    parfor i = 1:length(tr)
        estimates(i) = fminunc(@(v) -pop_motion_signal(v, spikes, cell_indices1, cell_indices2, cell_x_pos, tr(i), stop, run_opt.tau, run_opt.tol*.1), run_opt.trial_estimate_start, options);
        fprintf('for trial %d, the estimated speed was %d', i, estimates(i))
    end
    
    % stop parallel pool
    delete(poolobj);
    
    % save estimates
    save(sprintf('/home/vision/Dropbox/Lab/Development/matl ab-standard/private/malcolm/results/%s/downsample_%s_data_run_%02d_config_%d.mat', run_opt.data_set, run_opt.cell_type, run_opt.data_run, run_opt.config_num), 'estimates')
end

ElapsedTime=toc;