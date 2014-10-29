% Script to calculate motion signals across different cell populations
% This script calculates the motion signals for a data run displays several
% nice guis for visualization.  Can set lots of options for different
% behavior. 
%
% Marvin Thielk 2013
% mthielk@salk.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% vary the options set in run_opt to turn different behaviors on and off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load -> reloads data
% data_set -> sets the data set to do the analysis on.  built in support
%             for '2007-03-27-1' and '2007-08-24-4'
% remote -> if true uses data on the network, otherwise, local data
% data_run -> which data run to load and analyze
% cell_type -> string indicating which of the 4 cell populations to look at
% config_num -> which configuration of the stimuli to use (out of 4)
% raster -> creates a plot where you can vary the cell number and see a
%           raster for all the trials to demonstrate the reliability of
%           signals
% trial_raster -> creates a plot where you can vary the trial number and
%                 see all the spiking cells arranged by x position of their
%                 receptive field
% trial_raster_shift -> creates a plot where you can vary the trial number
%                       and manually shift the cells spike timing to see
%                       them line up.
% manual_speed_tuning -> creates a plot that demonstrates the algorithm and
%                        allows you see the gaussians line up as you vary
%                        the velocity.  This is for pairwise motion signal.
% velocity_lim -> max velocity that auto_speed_tuning automatically
%                 calculates
% auto_speed_tuning -> for a pair of neurons calculates the motion signal
%                      at different velocities and plots it.  Useful to
%                      show the pairwise shape of the speed tuning cuve
% tau -> tuning parameter that dictates the width of the gaussians used
% pop_speed_tuning -> plots the tuning curve summed over all possible pairs
%                     of neurons. Parallelized.
% tol -> determines the various tolerances for the integration and
%        optimizations used.
% savefig -> if true the script tries to save the figures that take longer
%            to produce in a figs/data_set folder.  The folder must already
%            exist or it will throw an error.
% trial_num -> Which trial to use for pop_speed_tuning.  
%                   0 < trial_num <= ~50
% trial_estimate -> creates a histogram of the fminunc determined velocity
%                   estimate of all the trials.  Parallelized.
% auto_set -> automatically sets a few of the parameters that I know work
%             best for different data runs.  CAUTION: this overwrites
%             values in run_opt.
% trial_estimate_start -> the velocity to use as the initial guess in
%                         fminunc to estimate the velocity.
% data_run_plots -> same as trial_estimate but does it for each of the four
%                   cell types
% cell_types -> the cell types to include when running data_run_plots 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

run_opt.load = true; % T/F
run_opt.data_set = '2007-03-27-1';
%run_opt.data_set = '2007-08-24-4';
run_opt.remote = true; % T/F
run_opt.data_run = 18; % 12-19
run_opt.cell_type = 'On parasol'; % on/off parasol, on/off midget
run_opt.config_num = 1; % 1-4
run_opt.raster = true; % T/F
run_opt.trial_raster = true; % T/F
run_opt.trial_raster_shift = true; % T/F
run_opt.manual_speed_tuning = true; % T/F
run_opt.velocity_lim = 150; % >0
run_opt.auto_speed_tuning = true; % T/F
run_opt.tau = .01; % tuning parameter
run_opt.pop_speed_tuning = true; % T/F
run_opt.tol = 1e-3;
run_opt.savefig = true; % T/F
run_opt.trial_num = 1; % > 0
run_opt.trial_estimate = false; % T/F
run_opt.auto_set = false; % T/F -- note: overwrites run_opt params
run_opt.trial_estimate_start = 120;
run_opt.data_run_plots = true; % T/F
run_opt.cell_types = {'Off midget', 'Off parasol', 'On midget', 'On parasol'};

% just some presets so I don't have to change every param each time I
% change the data_run
if run_opt.auto_set
    if strcmp(run_opt.data_set, '2007-03-27-1')
        run_opt.cell_types = {'Off midget', 'Off parasol', 'On midget', 'On parasol'};
        if run_opt.data_run == 12
            run_opt.velocity_lim = 150;
            run_opt.config_num = 1;
            run_opt.trial_estimate_start = 110;
            run_opt.tol = 1e-3;
        elseif run_opt.data_run == 13
            run_opt.velocity_lim = 300;
            run_opt.config_num = 1;
            run_opt.trial_estimate_start = 203;
            run_opt.tol = 1e-4;
        elseif run_opt.data_run == 14
            run_opt.velocity_lim = 50;
            run_opt.config_num = 1;
            run_opt.trial_estimate_start = 14.6;
            run_opt.tol = 1e-2;
        elseif run_opt.data_run == 15
            run_opt.velocity_lim = 150;
            run_opt.config_num = 3;
            run_opt.trial_estimate_start = 110;
            run_opt.tol = 1e-3;
        elseif run_opt.data_run == 16
            run_opt.velocity_lim = 300;
            run_opt.config_num = 1;
            run_opt.trial_estimate_start = 203;
            run_opt.tol = 1e-4;
        elseif run_opt.data_run == 17
            run_opt.velocity_lim = 50;
            run_opt.config_num = 2;
            run_opt.trial_estimate_start = 14.6;
            run_opt.tol = 1e-2;
        elseif run_opt.data_run == 18
            run_opt.velocity_lim = 150;
            run_opt.config_num = 1;
            run_opt.trial_estimate_start = 110;
            run_opt.tol = 1e-3;
        elseif run_opt.data_run == 19
            run_opt.velocity_lim = 300;
            run_opt.config_num = 1;
            run_opt.trial_estimate_start = 203;
            run_opt.tol = 1e-4;
        end
    elseif strcmp(run_opt.data_set, '2007-08-24-4')
        run_opt.cell_types = {'Off parasol', 'On midget', 'On parasol'};
        if run_opt.data_run == 2
            run_opt.velocity_lim = 150;
            run_opt.config_num = 2;
            run_opt.trial_estimate_start = 100;
            run_opt.tol = 1e-3;
        elseif run_opt.data_run == 3
            run_opt.velocity_lim = 300;
            run_opt.config_num = 3;
            run_opt.trial_estimate_start = 200;
            run_opt.tol = 1e-4;
        elseif run_opt.data_run == 4
            run_opt.velocity_lim = 300;
            run_opt.config_num = 1;
            run_opt.trial_estimate_start = 200;
            run_opt.tol = 1e-4;
        elseif run_opt.data_run == 5
            run_opt.velocity_lim = 300;
            run_opt.config_num = 2;
            run_opt.trial_estimate_start = 200;
            run_opt.tol = 1e-4;
        elseif run_opt.data_run == 6
            run_opt.velocity_lim = 150;
            run_opt.config_num = 3;
            run_opt.trial_estimate_start = 100;
            run_opt.tol = 1e-3;
        elseif run_opt.data_run == 7
            run_opt.velocity_lim = 150;
            run_opt.config_num = 1;
            run_opt.trial_estimate_start = 100;
            run_opt.tol = 1e-3;
        elseif run_opt.data_run == 8
            run_opt.velocity_lim = 50;
            run_opt.config_num = 1;
            run_opt.trial_estimate_start = 12.5;
            run_opt.tol = 1e-2;
        elseif run_opt.data_run == 9
            run_opt.velocity_lim = 50;
            run_opt.config_num = 1;
            run_opt.trial_estimate_start = 12.5;
            run_opt.tol = 1e-2;
        elseif run_opt.data_run == 10
            run_opt.velocity_lim = 150;
            run_opt.config_num = 2;
            run_opt.trial_estimate_start = 50;
            run_opt.tol = 1e-2;
        elseif run_opt.data_run == 11
            run_opt.velocity_lim = 150;
            run_opt.config_num = 1;
            run_opt.trial_estimate_start = 50;
            run_opt.tol = 1e-2;
        end
    end
end

% makes sure export_fig is on the path
if exist('export_fig', 'file') == 7
    addpath export_fig
end

if run_opt.load % load data fresh

    clear datarun tr

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

tic;

if run_opt.raster || run_opt.trial_raster || ...
        run_opt.trial_raster_shift || run_opt.manual_speed_tuning || ...
        run_opt.auto_speed_tuning || run_opt.pop_speed_tuning || ...
        run_opt.trial_estimate || run_opt.data_run_plots
    clf; set(gcf, 'color', 'white');
    
    cell_indices1=get_cell_indices(datarun{1},{run_opt.cell_type});
    cell_indices2=get_cell_indices(datarun{2},{run_opt.cell_type});
    
    cell_x_pos = cellfun( @(X) X.mean(1), datarun{1}.vision.sta_fits);
    [~, cell_sort_idx] = sort(cell_x_pos(cell_indices1));
    
    cell_indices1 = cell_indices1(cell_sort_idx);
    cell_indices2 = cell_indices2(cell_sort_idx);
    
    start = 0;
    stop = mean(datarun{2}.triggers(2:2:end) - datarun{2}.triggers(1:2:end));

    tr=datarun{2}.triggers(1:2:end); % triggers mark the beginning and end
    t=find(datarun{2}.stimulus.trial_list==run_opt.config_num);
    tr=tr(t);
end

if run_opt.raster %raster

    k=1; kmin=1; kmax=length(cell_indices2); hk=loop_slider(k,kmin,kmax);

    while k
        if ~ishandle(hk)
            break
        end
        k=round(get(hk,'Value')); 
 
        psth_raster(start,stop,datarun{2}.spikes{cell_indices2(k)}',tr);
        title(sprintf('%d %.2f', datarun{2}.cell_ids(cell_indices2(k)), datarun{1}.vision.sta_fits{cell_indices1(k)}.mean(1) ))

        uiwait;
    end
end

if run_opt.trial_raster
    k=1; kmin=1; kmax=length(tr); hk=loop_slider(k,kmin,kmax);
    while k
        if ~ishandle(hk)
            break
        end
        k = round(get(hk, 'Value'));
        
        shifted_trial_raster(start, stop, datarun{2}.spikes, cell_indices2, tr(k))
        title(sprintf('trial number %d', k))
        
        uiwait;
    end
end

if run_opt.trial_raster_shift
    k=1; kmin=1; kmax=length(tr); hk=loop_slider_n(k,kmin,kmax,1);
    j=run_opt.trial_estimate_start; jmin=-run_opt.velocity_lim; jmax=run_opt.velocity_lim; hj=loop_slider_n(j,jmin,jmax,2);
    while true
        if ~ishandle(hk)
            break
        end
        k = round(get(hk, 'Value'));
        j = round(get(hj, 'Value'));
        
        trigger = tr(k) + (cell_x_pos(cell_indices1) - cell_x_pos(cell_indices1(1))) ./ j;
        shifted_trial_raster(start, stop, datarun{2}.spikes, cell_indices2, trigger)
        title(sprintf('trial number %d with putative speed %d stixels/sec', k, j))
        
        uiwait;
    end
end

if run_opt.manual_speed_tuning
    k=1; kmin=1; kmax=length(tr); hk=loop_slider_n(k,kmin,kmax,1);
    jmin=1; jmax=length(cell_indices2); j=round(.15 * jmax); hj=loop_slider_n(j,jmin,jmax,2);
    hmin=1; hmax=length(cell_indices2); h=round(.9 * hmax); hh=loop_slider_n(h,hmin,hmax,3);
    v=run_opt.trial_estimate_start; vmin=-run_opt.velocity_lim; vmax=run_opt.velocity_lim; hv=loop_slider_n(v,vmin,vmax,4);
    while true
        if ~ishandle(hk)
            break
        end
        k = round(get(hk, 'Value'));
        j = round(get(hj, 'Value'));
        h = round(get(hh, 'Value'));
        v = get(hv, 'Value');
        
        spks_1 = datarun{2}.spikes{cell_indices2(j)};
        spks_2 = datarun{2}.spikes{cell_indices2(h)};
        dx = cell_x_pos(cell_indices1(h)) - cell_x_pos(cell_indices1(j));
        [sig_str, flt_rsp1, flt_rsp2, flt_rsp1_shifted, flt_rsp2_shifted, spks_1_shifted, spks_2_shifted] = motion_signal(v, spks_1, spks_2, dx, tr(k), stop, run_opt.tau);
        t = linspace(0, stop, 200);
        
        lims = [0 stop 0 4];
        
        subplot(3, 1, 1)
        plot(t, flt_rsp1(t), 'b-', t, flt_rsp2_shifted(t), 'g--', spks_2_shifted, ones(size(spks_2_shifted)), 'g.')
        title(sprintf('motion signal strength between neurons %d and %d in trial number %d', cell_indices2(j), cell_indices2(h), k))
        legend('filtered response of 1st cell', 'filtered response of 2nd cell shifted')
        axis(lims)
        
        subplot(3, 1, 2)
        plot(t, flt_rsp2(t), 'g-', t, flt_rsp1_shifted(t), 'b--', spks_1_shifted, ones(size(spks_1_shifted)), 'b.')
        title(sprintf('velocity = %d   and dx = %d', v, dx))
        legend('filtered response of 2nd cell', 'filtered response of 1st cell shifted')
        axis(lims)
        
        subplot(3, 1, 3)
        plot(t, flt_rsp1(t) .* flt_rsp2_shifted(t) - flt_rsp2(t) .* flt_rsp1_shifted(t))
        title(sprintf('signal strength: %d',sig_str))
        
        uiwait;
    end
end

if run_opt.auto_speed_tuning
    k=1; kmin=1; kmax=length(tr); hk=loop_slider_n(k,kmin,kmax,1);
    jmin=1; jmax=length(cell_indices2); j=round(.15 * jmax); hj=loop_slider_n(j,jmin,jmax,2);
    hmin=1; hmax=length(cell_indices2); h=round(.9 * hmax); hh=loop_slider_n(h,hmin,hmax,3);
    while true
        if ~ishandle(hk)
            break
        end
        k = round(get(hk, 'Value'));
        j = round(get(hj, 'Value'));
        h = round(get(hh, 'Value'));
        
        v = linspace(1, run_opt.velocity_lim);
        
        sig_str = zeros(size(v));
        for i = 1:length(v)
            spks_1 = datarun{2}.spikes{cell_indices2(j)};
            spks_2 = datarun{2}.spikes{cell_indices2(h)};
            dx = cell_x_pos(cell_indices1(h)) - cell_x_pos(cell_indices1(j));
            sig_str(i) = motion_signal(v(i), spks_1, spks_2, dx, tr(k), stop, run_opt.tau);
        end
        
        subplot(2,1,1)
        plot(v, sig_str)
        title(sprintf('motion signal strength between neurons %d and %d in trial number %d', cell_indices2(j), cell_indices2(h), k))
        xlabel('velocity')
        ylabel('net rightward motion signal')
        
        t = linspace(tr(k), tr(k)+stop, 500);
        trial_spks1 = spks_1(spks_1 >= tr(k) & spks_1 <= (tr(k) + stop));
        trial_spks2 = spks_2(spks_2 >= tr(k) & spks_2 <= (tr(k) + stop));
        flt_rsp1 = filtered_response(trial_spks1, run_opt.tau);
        flt_rsp2 = filtered_response(trial_spks2, run_opt.tau);
        subplot(4,1,3)
        plot(t, flt_rsp1(t), 'b', trial_spks1, ones(size(trial_spks1)), 'k.')
        lims = [tr(k) (tr(k) + stop) 0 max(max(flt_rsp1(t)), max(flt_rsp2(t))) * 1.1];
        axis(lims);
        title(sprintf('Average ISI: %d', (trial_spks1(end) - trial_spks1(1)) / (length(trial_spks1) - 1)));
        xlabel('time')
        
        subplot(4,1,4)
        plot(t, flt_rsp2(t), 'b', trial_spks2, ones(size(trial_spks2)), 'r.')
        axis(lims);
        title(sprintf('Average ISI: %d', (trial_spks2(end) - trial_spks2(1)) / (length(trial_spks2) - 1)));
        xlabel('time')
        
        uiwait;
    end
end

if run_opt.pop_speed_tuning
%     if matlabpool('size') <= 0
%         matlabpool
%     end
    
    v = linspace(1, run_opt.velocity_lim, 50);
    
    sig_str = zeros(size(v));
    for i = 1:length(v)
        sig_str(i) = pop_motion_signal(v(i), datarun{2}.spikes, cell_indices1, cell_indices2, cell_x_pos, tr(run_opt.trial_num), stop, run_opt.tau, run_opt.tol*.1);
        
        fprintf('*')
    end
    fprintf('\n')
    plot(v, sig_str)
    xlabel('velocity')
    ylabel('net rightward motion signal')
    title(sprintf('motion signal strength  in trial number %d', run_opt.trial_num))
%     if run_opt.savefig
%         export_fig(sprintf('figs/%s/%s_data_run_%02d_config_%d_trial_%d', run_opt.data_set, run_opt.cell_type, run_opt.data_run, run_opt.config_num, run_opt.trial_num), '-png', '-r300', '-painters')
%     end
end

if run_opt.trial_estimate
    if matlabpool('size') <= 0
        matlabpool
    end
    options = optimset('Display', 'iter', 'TolFun', run_opt.tol , 'MaxFunEvals', 30, 'LargeScale', 'off');
    estimates = zeros(size(tr));
    parfor i = 1:length(tr)
        estimates(i) = fminunc(@(v) -pop_motion_signal(v, datarun{2}.spikes, cell_indices1, cell_indices2, cell_x_pos, tr(i), stop, run_opt.tau, run_opt.tol*.1), run_opt.trial_estimate_start, options);
        fprintf('for trial %d, the estimated speed was %d', i, estimates(i))
    end
    figure()
    hist(estimates)
    xlabel('speed estimate (stixels/sec)')
    ylabel('trials')
    title(sprintf('%s data run %d config %d', run_opt.cell_type, run_opt.data_run, run_opt.config_num))
    if run_opt.savefig
        export_fig(sprintf('figs/%s/%s_data_run_%02d_config_%d', run_opt.data_set, run_opt.cell_type, run_opt.data_run, run_opt.config_num), '-png', '-r300', '-painters')
    end
    save(sprintf('data/%s/%s_data_run_%02d_config_%d.mat', run_opt.data_set, run_opt.cell_type, run_opt.data_run, run_opt.config_num), 'estimates')
end

if run_opt.data_run_plots
    if matlabpool('size') <= 0
        matlabpool
    end
    options = optimset('Display', 'iter', 'TolFun', run_opt.tol, 'MaxFunEvals', 30, 'LargeScale', 'off');
    for j=1:length(run_opt.cell_types)
        cell_type = run_opt.cell_types{j};
        cell_indices1=get_cell_indices(datarun{1},{cell_type});
        cell_indices2=get_cell_indices(datarun{2},{cell_type});
        
        cell_x_pos = cellfun( @(X) X.mean(1), datarun{1}.vision.sta_fits);
        [~, cell_sort_idx] = sort(cell_x_pos(cell_indices1));
        
        cell_indices1 = cell_indices1(cell_sort_idx);
        cell_indices2 = cell_indices2(cell_sort_idx);
        
        estimates = zeros(size(tr));
        parfor i = 1:length(tr)
            estimates(i) = fminunc(@(v) -pop_motion_signal(v, datarun{2}.spikes, cell_indices1, cell_indices2, cell_x_pos, tr(i), stop, run_opt.tau, run_opt.tol*.1), run_opt.trial_estimate_start, options);
            fprintf('for %s cells, in trial %d, the estimated speed was %d', cell_type, i, estimates(i))
        end
        figure()
        hist(estimates)
        xlabel('speed estimate (stixels/sec)')
        ylabel('trials')
        title(sprintf('%s data run %d config %d', cell_type, run_opt.data_run, run_opt.config_num))
        if run_opt.savefig
            export_fig(sprintf('figs/%s/%s_data_run_%02d_config_%d', run_opt.data_set, cell_type, run_opt.data_run, run_opt.config_num), '-png', '-r300', '-painters')
        end
        save(sprintf('data/%s/%s_data_run_%02d_config_%d.mat', run_opt.data_set, cell_type, run_opt.data_run, run_opt.config_num), 'estimates')
    end
end

toc

matlabpool close
