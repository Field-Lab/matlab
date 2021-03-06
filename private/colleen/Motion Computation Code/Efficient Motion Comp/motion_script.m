% This program will compute two types of rasters as well as the speed
% estimate from the population signal when two populations of cells are
% grouped together (like ON and OFF parasol cells)
% It can also downsample spikes but I haven't explored that part


% Saves a file containing the speed estimate for each trial to the specified
% folder
% DATA PARAMETERS


% This is a wrapper for motion_script_colleen that allows it to be called
% as a function.
% The parameters data_set, data_run, config_num, cell_type, vel as set in
% runMotionScriptColleen.

% Inputs:
% data_set: Date for run
% data_run: Run number
% config_num: Which stimulus condition to run
% cell_type: Which cell type to run
% vel: true speed

% Saves a file containing the speed estimate for each trial to the specified
% folder

function [estimates] = motion_script(varargin)
% USER SET PARAMETERS
% NUMERICAL PARAMETERS
run_opt.courseIter = 6;
run_opt.tau = .01; % tuning parameter
run_opt.tol = 1e-4;

run_opt.data_set_vec = {'2007-08-24-4'};
run_opt.data_run_vec =10; % 12-}19 for 2007-03-27, 2-11 for 2007-08-24, 13-17 for 2005-04-26
run_opt.direction_vec = {'right'}; % right or left
run_opt.contrast_vec = {'bright'};  %bright or dark
run_opt.velocity_exp_vec = 96; % >0
run_opt.cell_type_vec{1} = {'On parasol'}; % on/off parasol, on/off midget
%     run_opt.cell_type_vec{2} = {'On midget'}; % on/off parasol, on/off midget


if length(varargin) == 1
    run_opt.data_set_vec = {varargin{1}};
elseif length(varargin) == 2
    run_opt.data_set_vec = varargin{1};
    run_opt.data_run_vec = varargin{2}; % 12-19 for 2007-03-27, 2-11 for 2007-08-24, 13-17 for 2005-04-26
elseif length(varargin) == 3
    run_opt.data_set_vec = {varargin{1}};
    run_opt.data_run_vec = varargin{2}; % 12-19 for 2007-03-27, 2-11 for 2007-08-24, 13-17 for 2005-04-26
    run_opt.direction_vec = {varargin{3}};
elseif length(varargin) == 4
    run_opt.data_set_vec = {varargin{1}};
    run_opt.data_run_vec = varargin{2}; % 12-19 for 2007-03-27, 2-11 for 2007-08-24, 13-17 for 2005-04-26
    run_opt.direction_vec = {varargin{3}};
    run_opt.contrast_vec = {varargin{4}};
elseif length(varargin) == 5
    run_opt.data_set_vec = {varargin{1}};
    run_opt.data_run_vec = varargin{2}; % 12-19 for 2007-03-27, 2-11 for 2007-08-24, 13-17 for 2005-04-26
    run_opt.direction_vec = {varargin{3}};
    run_opt.contrast_vec = {varargin{4}};
    run_opt.velocity_exp_vec = varargin{5}; % >0
elseif length(varargin) == 6
    run_opt.data_set_vec = {varargin{1}};
    run_opt.data_run_vec = varargin{2}; % 12-19 for 2007-03-27, 2-11 for 2007-08-24, 13-17 for 2005-04-26
    run_opt.direction_vec = {varargin{3}};
    run_opt.contrast_vec = {varargin{4}};
    run_opt.velocity_exp_vec = varargin{5}; % >0
    run_opt.cell_type_vec{1} = {varargin{6}}; % on/off parasol, on/off midget
elseif length(varargin) == 7
    run_opt.data_set_vec = {varargin{1}};
    run_opt.data_run_vec = varargin{2}; % 12-19 for 2007-03-27, 2-11 for 2007-08-24, 13-17 for 2005-04-26
    run_opt.direction_vec = {varargin{3}};
    run_opt.contrast_vec = {varargin{4}};
    run_opt.velocity_exp_vec = varargin{5}; % >0
    run_opt.cell_type_vec{1} = {varargin{6}}; % on/off parasol, on/off midget
    run_opt.cell_type2_vec{2} = {varargin{7}}; % on/off parasol, on/off midget
end


for i = 1:length(run_opt.data_run_vec)
    
    run_opt.data_set = run_opt.data_set_vec{i};
    run_opt.data_run = run_opt.data_run_vec(i); % 12-19 for 2007-03-27, 2-11 for 2007-08-24, 13-17 for 2005-04-26
    run_opt.cell_type{1} = run_opt.cell_type_vec{i}; % on/off parasol, on/off midget
    run_opt.direction = run_opt.direction_vec{i}; % on/off parasol, on/off midget
    run_opt.contrast = run_opt.contrast_vec{i}; % on/off parasol, on/off midget
    if strcmp(run_opt.contrast, 'bright')
        run_opt.contrast = [0.48, 0.48, 0.48];
    else
        run_opt.contrast = [-0.48, -0.48, -0.48];
    end
    
    run_opt.velocity_exp = run_opt.velocity_exp_vec(i);
    if exist('run_opt.cell_type2_vec')
        run_opt.cell_type{2} = run_opt.cell_type2_vec{i}; % on/off parasol, on/off midget
    end
    
    % DATA PARAMETERS
    run_opt.cell_types = {'Off midget', 'Off parasol', 'On midget', 'On parasol'};
    run_opt.save = 0;
    
    % ANALYSES TO RUN
    run_opt.load = true; % T/F
    run_opt.auto_set = false; % T/F -- note: overwrites run_opt params
    run_opt.downsample_spikes = false; % must run on bertha
    run_opt.raster = false; % T/F
    run_opt.rasterPerTrial = true; % T/F
    run_opt.trial_estimate = true; % T/F
    
    
    if strcmp(run_opt.data_set, '2007-03-27-1')
        run_opt.stx = 8;
    elseif strcmp(run_opt.data_set, '2007-08-24-4')
        run_opt.stx = 10;
    else
        run_opt.stx = 10;
    end
    
    
    % Get computer identification
    sid = '';
    ni = java.net.NetworkInterface.getNetworkInterfaces;
    while ni.hasMoreElements
        addr = ni.nextElement.getHardwareAddress;
        if ~isempty(addr)
            addrStr = dec2hex(int16(addr)+128);
            sid = [sid, '.', reshape(addrStr,1,2*length(addr))];
        end
    end
    run_opt.sid = sid;
    % Bertha: '.163C7602AD5C';
    % My computer: '.CE033C0CF908'
    
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
    
    % Get the right configuration number
    run_opt.config_num = findStimType(datarun, run_opt.contrast, run_opt.direction);
    
    % Get indicies for each cell type
    for type  = 1:size(run_opt.cell_type,2)
        % Gets the indicies used by vision of the particular cell type
        if run_opt.raster || run_opt.trial_estimate || run_opt.rasterPerTrial
            
            
            % Get indices for specified cell type and order by RF position
            cell_indices1{type}=get_cell_indices(datarun{1},{run_opt.cell_type{type}});
            cell_indices2{type}=get_cell_indices(datarun{2},{run_opt.cell_type{type}});
            cell_x_pos = cellfun( @(X) X.mean(1), datarun{1}.vision.sta_fits); % x axis position of all STA cells
            [~, cell_sort_idx{type}] = sort(cell_x_pos(cell_indices1{type})); % indicies of how to sort
            
            %cell_indices sorted by their x coordinate of the RF from the STA
            cell_indices1{type} = cell_indices1{type}(cell_sort_idx{type}); % cell_indices1 is now indexes in order from lowest to highest firing rate
            cell_indices2{type} = cell_indices2{type}(cell_sort_idx{type});
            
            % Find trial start and stop times
            start = 0;
            stop = mean(datarun{2}.triggers(2:2:end) - datarun{2}.triggers(1:2:end));
            tr=datarun{2}.triggers(1:2:end); % all start triggers
            t=find(datarun{2}.stimulus.trial_list==run_opt.config_num); %find the times when all the stimulus type 2 starts
            tr=tr(t);
        end
    end
    % Grouped pooled together
    if size(run_opt.cell_type,2) == 1
        cell_indices1 = [cell_indices1{1}];
        cell_indices2 = [cell_indices2{1}];
        %     cell_x_pos = [cell_x_pos{1}];
        [~, cell_sort_idx] = sort(cell_x_pos(cell_indices1)); % indicies of how to sort
        
    else
        
        cell_indices1 = [cell_indices1{1}, cell_indices1{2}];
        cell_indices2 = [cell_indices2{1}, cell_indices2{2}];
        %     cell_x_pos = [cell_x_pos{1}, cell_x_pos{2}];
        [~, cell_sort_idx] = sort(cell_x_pos(cell_indices1)); % indicies of how to sort
    end
    
    %cell_indices sorted by their x coordinate of the RF from the STA
    cell_indices1= cell_indices1(cell_sort_idx); % cell_indices1 is now indexes in order from lowest to highest firing rate
    cell_indices2 = cell_indices2(cell_sort_idx);
    
    
    
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
                if strcmp(run_opt.cell_type{type}, 'On parasol')
                    for i=cell_indices2
                        % fractional downsampling:
                        datarun{2}.spikes{i}=datasample(datarun{2}.spikes{i},round(length(datarun{2}.spikes{i})*n2/n1),'Replace',false);
                    end
                else
                    run_opt.trial_estimate = false; % if no downsampling, don't calculate estimates
                end
            elseif n2>n1
                if strcmp(run_opt.cell_type{type}, 'On midget')
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
                if strcmp(run_opt.cell_type{type}, 'On parasol')
                    for i=cell_indices2
                        % fractional downsampling
                        datarun{2}.spikes{i}=datasample(datarun{2}.spikes{i},round(length(datarun{2}.spikes{i})*n2/n1),'Replace',false);
                    end
                else
                    run_opt.trial_estimate = false; % if no downsampling, don't calculate estimates
                end
            elseif n2>n1
                if strcmp(run_opt.cell_type{type}, 'On midget')
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
    
    
    
    % Plot one cell on all trials
    if run_opt.raster %raster
        cell_indices1 = [cell_indices1{1}, cell_indices1{2}];
        cell_indices2 = [cell_indices2{1}, cell_indices2{2}];
        cell_x_pos = [cell_x_pos{1}, cell_x_pos{2}];
        [~, cell_sort_idx] = sort(cell_x_pos(cell_indices1)); % indicies of how to sort
        
        %cell_indices sorted by their x coordinate of the RF from the STA
        cell_indices1= cell_indices1(cell_sort_idx); % cell_indices1 is now indexes in order from lowest to highest firing rate
        cell_indices2 = cell_indices2(cell_sort_idx);
        k=1; kmin=1; kmax=length(cell_indices2); hk=loop_slider(k,kmin,kmax);
        
        while k
            if ~ishandle(hk)
                break % script breaks until figure is closed
            end
            k=round(get(hk,'Value'));
            % Takes in start and stop time (0-0.7274)
            % Spikes of the cell with the lowest firing rate first
            % start time of each stimulus type trigger
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
    
    % Plot all cells on one trial
    % ONLY WORKS FOR RIGHT MOVING BAR
    if run_opt.rasterPerTrial
        
        toPlot = cell(1,length(cell_indices2));
        
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
         
        for counter = 1:length(cell_indices2)
            psth_r = psth_raster_noPlotting(start,stop,datarun{2}.spikes{cell_indices2(counter)}',tr, '','',counter);
            posThisCell = datarun{1}.vision.sta_fits{cell_indices1(counter)}.mean(1);
            
            posFarthestCell = datarun{1}.vision.sta_fits{cell_indices1(1)}.mean(1);
            
            
            cellNumber = datarun{2}.cell_ids(cell_indices2(counter));
            % Title is the cell id according to vision and the mean firing rate
            %          [psth, bins] = get_psth(datarun{2}.spikes{cell_indices2(counter)}, tr, 'plot_hist', true)
        
%             for trialNum = 1:length(t)
%                 [x,y] = find(psth_r == trialNum-1);
                toPlot{counter}= [[psth_r(:,1), repmat(cellNumber, size(psth_r,1),1), repmat(posThisCell, size(psth_r,1),1)]];
             hold on
                plot(toPlot{counter}(:,1),toPlot{counter}(:,3),Color);
      
        end
        
        k=1; kmin=1; kmax=length(t); hk=loop_slider(k,kmin,kmax);
        
        %         while k
        if ~ishandle(hk)
            break % script breaks until figure is closed
        end
        k=round(get(hk,'Value'));
        y_scale = 1;
       
        %             title({run_opt.cell_type{type}, [run_opt.data_set, ' Run ', num2str(run_opt.data_run)],'Bright bar moving right', sprintf(' Trial Number %d',  k)})
        xlabel('time (ms)');
        ylabel('Cell''s centroid distance from reference');
        %             uiwait;
    end
    
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
    k=1; kmin=1; kmax=length(t); hk=loop_slider(k,kmin,kmax);
    
    
    k=round(get(hk,'Value'));
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

% send email when done
% gmail('crhoades227@gmail.com', sprintf('Done with %s %s_data_run_%02d_config_%d_darkright_newmethod',run_opt.data_set, run_opt.cell_type, run_opt.data_run, run_opt.config_num))


end
