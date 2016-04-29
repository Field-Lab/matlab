%% Process moving bar; photons
clear

% STA PATH
run_opt{1}.data_set = '2016-04-21-1';
run_opt{1}.data_run = 'data000';


% MOVING BAR PATH
run_opt{2}.data_set = '2016-04-21-1';
run_opt{2}.data_run = 'data009-from-data000';

location = 'Data';

% Where to save the data
filepath= ['/Volumes/Lab/Users/crhoades/Moving_Bar/', run_opt{2}.data_set, '/', run_opt{2}.data_run, '/'];

%%%%%%%%%%%%%%%%%%%%%%%%%% END INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% STA PATH
datarun{1}.names.rrs_params_path=['/Volumes/Analysis/', run_opt{1}.data_set, '/', run_opt{1}.data_run,'/', run_opt{1}.data_run,  '.params'];
datarun{1}.names.rrs_params_path=['/Volumes/Analysis/', run_opt{1}.data_set, '/', run_opt{1}.data_run, '/', run_opt{1}.data_run, '.params'];
datarun{1}.names.rrs_sta_path=['/Volumes/Analysis/', run_opt{1}.data_set, '/', run_opt{1}.data_run, '/', run_opt{1}.data_run, '.sta'];


%% MOVING BAR PATH
datarun{2}.names.rrs_neurons_path=['/Volumes/Analysis/', run_opt{2}.data_set, '/', run_opt{2}.data_run, '/', run_opt{2}.data_run, '.neurons'];
datarun{2}.names.stimulus_path=['/Volumes/', location, '/', run_opt{2}.data_set, '/Visual/s_files/','s', run_opt{2}.data_run(6:7),'.txt'];
datarun{2}.names.rrs_params_path=['/Volumes/Analysis/', run_opt{2}.data_set, '/', run_opt{2}.data_run, '/', run_opt{2}.data_run, '.params'];
%
% datarun.names.rrs_neurons_path=['/Volumes/Analysis/', run_opt.data_set, '/', run_opt.data_run, '/', run_opt.data_run, '.neurons'];
% datarun.names.stimulus_path=['/Volumes/',location,'/', run_opt.data_set, '/Visual/s', run_opt.data_run(end-1:end)];
opt=struct('verbose',1,'load_sta', 1,'load_params',1,'load_neurons',1,'load_obvius_sta_fits',true);

datarun{1}=load_data(datarun{1},opt);
datarun{2}=load_data(datarun{2},opt);

% datarun=map_cell_types(datarun, struct('map',[1 2],'verbose',true));

% If the trigger interval in LabView is set long enought (~6 seconds for 5
% second stimuli), then this trigger_iti_thr should be fine.
datarun{2}=load_stim(datarun{2},'correction_incomplet_run', 0, 'trigger_iti_thr', 0.01); % manually set threshold until got right number of triggers



% for i = 1:length(run_opt.data_run_vec)

%     run_opt.data_set = run_opt.data_set_vec{i};
%     run_opt.data_run = run_opt.data_run_vec(i); % 12-19 for 2007-03-27, 2-11 for 2007-08-24, 13-17 for 2005-04-26
% run_opt{2}.cell_type{1} = 'ON parasol'; % on/off parasol, on/off midget
run_opt{2}.direction = 'left'; % on/off parasol, on/off midget
run_opt{2}.contrast = 'dark'; % on/off parasol, on/off midget
if strcmp(run_opt{2}.contrast, 'bright')
    run_opt{2}.contrast = [0.48, 0.48, 0.48];
else
    run_opt{2}.contrast = [-0.48, -0.48, -0.48];
end

run_opt{2}.delta = 8;
% if exist('run_opt.cell_type2_vec')
%     run_opt{2}.cell_type{2} = run_opt{2}.cell_type2_vec{i}; % on/off parasol, on/off midget
% end

% DATA PARAMETERS
run_opt{2}.cell_types = {'ON large'};
run_opt{2}.save = 0;

% ANALYSES TO RUN
run_opt{2}.load = true; % T/F
run_opt{2}.auto_set = false; % T/F -- note: overwrites run_opt params
run_opt{2}.downsample_spikes = false; % must run on bertha
run_opt{2}.raster = true; % T/F
run_opt{2}.rasterPerTrial = true; % T/F
run_opt{2}.trial_estimate = false; % T/F


run_opt{2}.stx = 8;


% Get computer identification
%     sid = '';
%     ni = java.net.NetworkInterface.getNetworkInterfaces;
%     while ni.hasMoreElements
%         addr = ni.nextElement.getHardwareAddress;
%         if ~isempty(addr)
%             addrStr = dec2hex(int16(addr)+128);
%             sid = [sid, '.', reshape(addrStr,1,2*length(addr))];
%         end
%     end
%     run_opt.sid = sid;
% Bertha: '.163C7602AD5C';
% My computer: '.CE033C0CF908'

tic;

% Auto set parameters if flag set to true
%     if run_opt.auto_set
%         [run_opt.cell_types, run_opt.velocity_lim, run_opt.config_num, run_opt.trial_estimate_start, run_opt.tol] =...
%             auto_set_params(run_opt.data_set, run_opt.data_run);
%     end

% Load data
%     if run_opt.load
%         clear datarun tr
%         % datarun{1} has vision info (sta fits)
%         % datarun{2} has cell_ids, spikes, triggers
%         if strcmp(run_opt.data_set, '2007-03-27-1')
%             datarun{1}.names.rrs_params_path='/Volumes/Analysis/2007-03-27-1/data011-nwpca/data011-nwpca.params';
%             datarun{2}.names.rrs_neurons_path=sprintf('/Volumes/Analysis/2007-03-27-1/data%03d-from-data011-nwpca/data%03d-from-data011-nwpca.neurons', run_opt.data_run, run_opt.data_run);
%             datarun{2}.names.stimulus_path=sprintf('/Volumes/Analysis/2007-03-27-1/stimuli/s%02d', run_opt.data_run);
%         elseif strcmp(run_opt.data_set, '2007-08-24-4')
%             datarun{1}.names.rrs_params_path='/Volumes/Analysis/2007-08-24-4/data001-nwpca/data001-nwpca.params';
%             datarun{2}.names.rrs_neurons_path=sprintf('/Volumes/Analysis/2007-08-24-4/data%03d-from-data001-nwpca/data%03d-from-data001-nwpca.neurons', run_opt.data_run, run_opt.data_run);
%             datarun{2}.names.stimulus_path=sprintf('/Volumes/Analysis/2007-08-24-4/Stimuli/s%02d', run_opt.data_run);
%         elseif strcmp(run_opt.data_set, '2005-04-26-0')
%             datarun{1}.names.rrs_params_path='/Volumes/Analysis/2005-04-26-0/data';
%         end
%         opt=struct('verbose',1,'load_params',1,'load_neurons',1,'load_obvius_sta_fits',true);
%         datarun=load_data(datarun,opt);
%         datarun=map_cell_types(datarun, struct('map',[1 2],'verbose',true));
%         datarun{2}=load_stim(datarun{2},'correction_incomplet_run', 0);
%     end

% Get the right configuration number
%     run_opt.config_num = findStimType(datarun, run_opt.delta, run_opt.direction);



%parameters are direction and delta
directions = [];
deltas = [];
for i = 1:size(datarun{2}.stimulus.trials,2)
    if sum(datarun{2}.stimulus.trials(i).direction == directions) == 0
        directions = [directions; datarun{2}.stimulus.trials(i).direction];
    end
    
    if sum(datarun{2}.stimulus.trials(i).delta == deltas) == 0
        deltas= [deltas; datarun{2}.stimulus.trials(i).delta];
    end
    
end

direction = directions(1);
delta = deltas(3);

for i = 1:size(datarun{2}.stimulus.combinations,2)
    if (datarun{2}.stimulus.combinations(i).direction == direction && datarun{2}.stimulus.combinations(i).delta == delta)
        trial_type = i;
        break;
    end
end

% Get indicies for each cell type
for type  = 1:size(run_opt{2}.cell_types,2)
    % Gets the indicies used by vision of the particular cell type
    if run_opt{2}.raster || run_opt{2}.trial_estimate || run_opt{2}.rasterPerTrial
        
        
        % Get indices for specified cell type and order by RF position
        cell_indices=get_cell_indices(datarun{1},run_opt{2}.cell_types(type));
        cell_x_pos = cellfun( @(X) X.mean(1), datarun{1}.vision.sta_fits); % x axis position of all STA cells
        [~, cell_sort_idx] = sort(cell_x_pos(cell_indices)); % indicies of how to sort
        
        %cell_indices sorted by their x coordinate of the RF from the STA
        cell_indices_sorted = cell_indices(cell_sort_idx); % cell_indices1 is now indexes in order from lowest to highest firing rate
        
        % Find trial start and stop times
        start = 0;
        stop = mean(datarun{2}.triggers(2:2:end) - datarun{2}.triggers(1:2:end));
        tr=datarun{2}.triggers(1:2:end); % all start triggers
        t=find(datarun{2}.stimulus.trial_list==trial_type); %find the times when all the stimulus type 2 starts
        tr=tr(t);
    end
end
% Grouped pooled together
% if size(run_opt{2}.cell_type,2) == 1
%     cell_indices1 = [cell_indices1{1}];
%     %     cell_x_pos = [cell_x_pos{1}];
%     [~, cell_sort_idx] = sort(cell_x_pos(cell_indices1)); % indicies of how to sort
%
% else
%
%     cell_indices1 = [cell_indices1{1}, cell_indices1{2}];
%     %     cell_x_pos = [cell_x_pos{1}, cell_x_pos{2}];
%     [~, cell_sort_idx] = sort(cell_x_pos(cell_indices1)); % indicies of how to sort
% % end
%
% %cell_indices sorted by their x coordinate of the RF from the STA
% cell_indices1= cell_indices1(cell_sort_idx); % cell_indices1 is now indexes in order from lowest to highest firing rate


%
% % downsample spikes
% if run_opt{2}.downsample_spikes
%
%     onp_indices = get_cell_indices(datarun{2}, {'On parasol'});
%     numspikes1 = zeros(length(onp_indices),1);
%     for i=1:length(onp_indices)
%         numspikes1(i) = length(datarun{2}.spikes{onp_indices(i)});
%     end
%
%     onm_indices = get_cell_indices(datarun{2}, {'On midget'});
%     numspikes2 = zeros(length(onm_indices),1);
%     for i=1:length(onm_indices)
%         numspikes2(i) = length(datarun{2}.spikes{onm_indices(i)});
%     end
%
%     n1 = mean(numspikes1);
%     n2 = mean(numspikes2);
%
%     if strcmp(run_opt.data_set, '2007-03-27-1')
%         if n1>n2
%             if strcmp(run_opt.cell_type{type}, 'On parasol')
%                 for i=cell_indices2
%                     % fractional downsampling:
%                     datarun{2}.spikes{i}=datasample(datarun{2}.spikes{i},round(length(datarun{2}.spikes{i})*n2/n1),'Replace',false);
%                 end
%             else
%                 run_opt.trial_estimate = false; % if no downsampling, don't calculate estimates
%             end
%         elseif n2>n1
%             if strcmp(run_opt.cell_type{type}, 'On midget')
%                 for i=cell_indices2
%                     % fractional downsampling
%                     datarun{2}.spikes{i}=datasample(datarun{2}.spikes{i},round(length(datarun{2}.spikes{i})*n1/n2),'Replace',false);
%                 end
%             else
%                 run_opt.trial_estimate = false; % if no downsampling, don't calculate estimates
%             end
%         end
%     elseif strcmp(run_opt.data_set, '2007-08-24-4')
%         if n1>n2
%             if strcmp(run_opt.cell_type{type}, 'On parasol')
%                 for i=cell_indices2
%                     % fractional downsampling
%                     datarun{2}.spikes{i}=datasample(datarun{2}.spikes{i},round(length(datarun{2}.spikes{i})*n2/n1),'Replace',false);
%                 end
%             else
%                 run_opt.trial_estimate = false; % if no downsampling, don't calculate estimates
%             end
%         elseif n2>n1
%             if strcmp(run_opt.cell_type{type}, 'On midget')
%                 for i=cell_indices2
%                     % fractional downsampling
%                     datarun{2}.spikes{i}=datasample(datarun{2}.spikes{i},round(length(datarun{2}.spikes{i})*n1/n2),'Replace',false);
%                 end
%             else
%                 run_opt.trial_estimate = false; % if no downsampling, don't calculate estimates
%             end
%         end
%     end
% end
% cell_indices = get_cell_indices(datarun{2}, run_opt{2}.cell_types{1})

% Plot one cell on all trials
if run_opt{2}.raster %raster
    %     cell_indices1 = [cell_indices1{1}, cell_indices1{2}];
    %     cell_x_pos = [cell_x_pos{1}, cell_x_pos{2}];
    %     [~, cell_sort_idx] = sort(cell_x_pos(cell_indices1)); % indicies of how to sort
    
    %cell_indices sorted by their x coordinate of the RF from the STA
    for k =1:length(cell_indices)
        
        %     cell_indices1= cell_indices1(cell_sort_idx{k}); % cell_indices1 is now indexes in order from lowest to highest firing rate
        
        % Takes in start and stop time (0-0.7274)
        % Spikes of the cell with the lowest firing rate first
        % start time of each stimulus type trigger
        % Finds the spikes that happened on a cell from stimulus onset to end
        % Plot those spike times on the x axis versus the trial number on the y
        % axis
        % If tracking motion, the cell should respond to the bar at the same
        % time on every trial
        psth_r =  psth_raster(start,stop,datarun{2}.spikes{cell_indices_sorted(k)}',tr);
%                 figure;
%                     rasterplot(psth_r(:,1),max(psth_r(:,2)+1) ,3)
        [psth, bins] = get_psth(datarun{2}.spikes{cell_indices_sorted(k)}', tr, 'stop', stop, 'start', start, 'bin_size', 0.01);
%         figure;
%         plot(psth_r(:,1),psth_r(:,2),'k.');%, 'MarkerSize',10
%         axis([start*1000 stop*1000 0 length(tr)]);
%         % Title is the cell id according to vision and the mean firing rate
        title(sprintf('%d %.2f', datarun{2}.cell_ids(cell_indices(k)), datarun{1}.vision.sta_fits{cell_indices(k)}.mean(1) ))
    end
    
end

% Plot all cells on one trial
% ONLY WORKS FOR RIGHT MOVING BAR
if run_opt.rasterPerTrial
    
    toPlot = cell(1,length(cell_indices));
    
    % Takes in start and stop time (0-0.7274)
    % Spikes of the cell with the lowest firing rate first
    % start time of each stimulus type 2 trigger
    % Finds the spikes that happened on a cell from stimulus onset to end
    % Plot those spike times on the x axis versus the trial number on the y
    % axis
    % If tracking motion, the cell should respond to the bar at the same
    % time on every trial
    psth_r =[];
    Color = ['k', '.'];
    figure
    hold on
    
    for counter = 1:length(cell_indices)
    
        psth_r = psth_raster_noPlotting(start,stop,datarun{2}.spikes{cell_indices(counter)}',tr, '','',counter);
        posThisCell = datarun{1}.vision.sta_fits{cell_indices(counter)}.mean(1);
        
        posFarthestCell = datarun{1}.vision.sta_fits{cell_indices(1)}.mean(1);
        
        
        cellNumber = datarun{2}.cell_ids(cell_indices(counter));
        % Title is the cell id according to vision and the mean firing rate
        %          [psth, bins] = get_psth(datarun{2}.spikes{cell_indices2(counter)}, tr, 'plot_hist', true)
        
        %             for trialNum = 1:length(t)
        %                 [x,y] = find(psth_r == trialNum-1);
        toPlot{counter}= [[psth_r(:,1), repmat(cellNumber, size(psth_r,1),1), repmat(posThisCell, size(psth_r,1),1)]];
        hold on
        plot(toPlot{counter}(:,1),toPlot{counter}(:,3),Color);
        
    end
    
    %     k=1; kmin=1; kmax=length(t); hk=loop_slider(k,kmin,kmax);
    
    %         while k
    %     if ~ishandle(hk)
    %         break % script breaks until figure is closed
    %     end
    %     k=round(get(hk,'Value'));
    %     y_scale = 1;
    %
    %             title({run_opt.cell_type{type}, [run_opt.data_set, ' Run ', num2str(run_opt.data_run)],'Bright bar moving right', sprintf(' Trial Number %d',  k)})
    %     xlabel('time (ms)');
    %     ylabel('Cell''s centroid distance from reference');
    %             uiwait;
end





if run_opt.trial_estimate
    
    toPlot = cell(1,length(t));
    
    % Takes in start and stop time (0-0.7274)
    % Spikes of the cell with the lowest firing rate first
    % start time of each stimulus type 2 trigger
    % Finds the spikes that happened on a cell from stimulus onset to end
    % Plot those spike times on the x axis versus the trial number on the y
    % axis
    % If tracking motion, the cell should respond to the bar at the same
    % time on every trial
    %         for counter = 1:length(cell_indices2)
    %             psth_r = psth_raster_noPlotting(start,stop,datarun{2}.spikes{cell_indices2(counter)}',tr(1));
    %             posThisCell = datarun{1}.vision.sta_fits{cell_indices1(counter)}.mean(1);
    %
    %             posFarthestCell = datarun{1}.vision.sta_fits{cell_indices1(1)}.mean(1);
    %
    %
    %             cellNumber = datarun{2}.cell_ids(cell_indices2(counter));
    %             % Title is the cell id according to vision and the mean firing rate
    %             %          [psth, bins] = get_psth(datarun{2}.spikes{cell_indices2(counter)}, tr, 'plot_hist', true)
    %             for trialNum = 1%:length(t)
    %                 [x,y] = find(psth_r == trialNum-1);
    %                 toPlot{trialNum}= [toPlot{trialNum}; [psth_r(x,1), repmat(cellNumber, length(x),1), repmat(posThisCell, length(x),1)]];
    %             end
    %         end
    %
    %     k=1; kmin=1; kmax=length(t); hk=loop_slider(k,kmin,kmax);
    %
    %
    %     k=round(get(hk,'Value'));
    y_scale = 1;
    Color = ['k', '.'];
    
    
    options = optimset('Display', 'iter', 'TolFun', run_opt.tol , 'MaxFunEvals', 60, 'LargeScale', 'off');
    spikes = datarun{2}.spikes;
    
    %Prior is +/-25% of expected value
    velocity = 186:0.5:198;%linspace(0.75*run_opt.velocity_exp, 1.25*run_opt.velocity_exp, run_opt.courseIter);
    strsig1 = zeros(1,length(velocity));
    
    % Run coarse error function to initialize velocity
    %         for i =1:length(tr)
    %             parfor j = 1:length(velocity)
    %                 v = velocity(j);
    %                 [strsig1(j)] = -pop_motion_signal(v, spikes, cell_indices1, cell_indices2, cell_x_pos, tr(i), stop, run_opt.tau, run_opt.tol, datarun, run_opt.direction, run_opt.stx);
    %             end
    %                   figure; plot(velocity, strsig1)
    %             i
    %             [x1,y1] = min(strsig1);
    %
    %             Initialize minimization
    %         end
    %
    % Find speed estimate
    end_trial = stop;
    
    for j =0:29%:length(tr)
        %         for
        vel = 1%:length(velocity)
        run_opt.trial_estimate_start(vel) = velocity(vel);
        
        start= 0;%j*end_trial/10;
        stop = (j+1)*end_trial/30;
        if j > 0
            [estimates(j+1)] = fminunc(@(v) -pop_motion_signal(v, spikes, cell_indices1, cell_indices2, cell_x_pos, tr(1), start,stop, run_opt.tau, run_opt.tol, datarun, run_opt.direction, run_opt.stx), estimates(j), options);
        else
            [estimates(j+1)] = fminunc(@(v) -pop_motion_signal(v, spikes, cell_indices1, cell_indices2, cell_x_pos, tr(1), start,stop, run_opt.tau, run_opt.tol, datarun, run_opt.direction, run_opt.stx), run_opt.velocity_exp_vec, options);
        end
        
        
        %             [trash(vel), spksShiftedRight(:, vel)] = pop_motion_signal_getSpikes(estimates(vel), spikes, cell_indices1, cell_indices2, cell_x_pos, tr(1), start,stop, run_opt.tau, run_opt.tol, datarun, run_opt.direction, run_opt.stx);
        
        
        
        
        %             estimates(vel) = -estimates(vel);
        
        %             fprintf('for trial %d, the estimated speed was %d', j, estimates(j))
        %         end
        
        
        
        %         figure; plot(velocity, estimates);
        %         [x,y] = min(estimates);
        %         vel_tot(j+1) = estimates(j+1)%velocity(y);
        %         spksShiftedRight = spksShiftedRight(:, y);
        %         for counter = 1:length(cell_indices2)
        %             posThisCell = datarun{1}.vision.sta_fits{cell_indices1(counter)}.mean(1);
        %
        %
        %             cellNumber = datarun{2}.cell_ids(cell_indices2(counter));
        %             % Title is the cell id according to vision and the mean firing rate
        %             %          [psth, bins] = get_psth(datarun{2}.spikes{cell_indices2(counter)}, tr, 'plot_hist', true)
        %             for trialNum = 1%:length(t)
        %
        %                 toPlot{trialNum}= [toPlot{trialNum}; [spksShiftedRight{counter}*1000, repmat(cellNumber, length(spksShiftedRight{counter}),1), repmat(posThisCell, length(spksShiftedRight{counter}),1)]];
        %             end
        %         end
        
        
        %         window(j+1) = estimates(vel)% velocity(y);
        %         figure;
        %         plot(toPlot{1}(:,1),toPlot{1}(:,3)*y_scale,Color);
        xlabel('time (ms)');
        ylabel('Cell''s centroid distance from reference');
        
        
    end
    x_axis = linspace(start,stop, length(estimates));
    figure; plot(x_axis, estimates', 'o-')
    title('Speed Estimate Given Spikes Up To Certain Time')
    xlabel('Time')
    ylabel('Speed estimate')
    
    run_opt.elapsedTime=toc;
    %save results
    %     foldername = sprintf('/Users/vision/Desktop/GitHub code repository/private/colleen/Results/%s/BrightRight/OnMandOnP', run_opt.data_set);
    %     filename = sprintf('/Users/vision/Desktop/GitHub code repository/private/colleen/Results/%s/BrightRight/OnMAndOnP/%s_data_run_%02d_config_%d.mat', run_opt.data_set, run_opt.cell_type{type}, run_opt.data_run, run_opt.config_num);
    if run_opt.save
        foldername = sprintf('/home/vision/Colleen/matlab/private/colleen/results/%s/DarkLeft/OnMandOnP', run_opt.data_set);
        filename = sprintf('/home/vision/Colleen/matlab/private/colleen/results/%s/DarkLeft/OnMandOnP/%s_data_run_%02d_config_%d.mat', run_opt.data_set, run_opt.cell_type{type}, run_opt.data_run, run_opt.config_num);
        
        if exist(foldername)
            save(filename, 'estimates', 'run_opt');
        else
            mkdir(foldername);
            save(filename, 'estimates', 'run_opt');
        end
    end
    
end

