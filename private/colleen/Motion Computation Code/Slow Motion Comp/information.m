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
% run_opt.data_set = '2007-03-27-1';
run_opt.data_set = '2007-08-24-4';
run_opt.remote = true; % T/F
run_opt.data_run = 3; % 12-19
run_opt.cell_type = 'On parasol'; % on/off parasol, on/off midget
run_opt.config_num = 1; % 1-4
run_opt.raster = false; % T/F
run_opt.trial_raster = false; % T/F
run_opt.trial_raster_shift = false; % T/F
run_opt.manual_speed_tuning = false; % T/F
run_opt.velocity_lim = 150; % >0
run_opt.auto_speed_tuning = false; % T/F
run_opt.tau = .01; % tuning parameter
run_opt.pop_speed_tuning = false; % T/F
run_opt.tol = 1e-3;
run_opt.savefig = true; % T/F
run_opt.trial_num = 1; % > 0
run_opt.trial_estimate = false; % T/F
run_opt.auto_set = false; % T/F -- note: overwrites run_opt params
run_opt.trial_estimate_start = 120;
run_opt.data_run_plots = false; % T/F
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
        if run_opt.remote 
            datarun{1}.names.rrs_params_path='/snle/analysis/2007-03-27-1/data011-nwpca/data011-nwpca.params';
            datarun{2}.names.rrs_neurons_path=sprintf('/snle/analysis/2007-03-27-1/data%03d-from-data011-nwpca/data%03d-from-data011-nwpca.neurons', run_opt.data_run, run_opt.data_run);
            datarun{2}.names.stimulus_path=sprintf('/braid/snle/analysis-archive/Experiments/Array/Analysis/2007-03-27-1/stimuli/s%02d', run_opt.data_run);
        else
            datarun{1}.names.rrs_params_path='/Data/2007-03-27-1/data011-nwpca/data011-nwpca.params';
            datarun{2}.names.rrs_neurons_path=sprintf('/Data/2007-03-27-1/data%03d-from-data011-nwpca/data%03d-from-data011-nwpca.neurons', run_opt.data_run, run_opt.data_run);
            datarun{2}.names.stimulus_path=sprintf('/braid/snle/analysis-archive/Experiments/Array/Analysis/2007-03-27-1/stimuli/s%02d', run_opt.data_run);
        end
    elseif strcmp(run_opt.data_set, '2007-08-24-4')
        if run_opt.remote 
            datarun{1}.names.rrs_params_path='/snle/analysis/2007-08-24-4/data001-nwpca/data001-nwpca.params';
            datarun{2}.names.rrs_neurons_path=sprintf('/snle/analysis/2007-08-24-4/data%03d-from-data001-nwpca/data%03d-from-data001-nwpca.neurons', run_opt.data_run, run_opt.data_run);
            datarun{2}.names.stimulus_path=sprintf('/braid/snle/analysis-archive/Experiments/Array/Analysis/2007-08-24-4/Stimuli/s%02d', run_opt.data_run);
        else
            datarun{1}.names.rrs_params_path='/Data/2007-08-24-4/data001-nwpca/data001-nwpca.params';
            datarun{2}.names.rrs_neurons_path=sprintf('/Data/2007-08-24-4/data%03d-from-data001-nwpca/data%03d-from-data001-nwpca.neurons', run_opt.data_run, run_opt.data_run);
            datarun{2}.names.stimulus_path=sprintf('/braid/snle/analysis-archive/Experiments/Array/Analysis/2007-08-24-4/Stimuli/s%02d', run_opt.data_run);
        end
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

cell_idxs = cell_indices2;
dt = .1;
bins = start:dt:stop;

bin_spikes = zeros(length(cell_idxs), length(bins)-1,length(tr));

for k = 1:length(tr)
spikes = datarun{2}.spikes
trigger = tr(k);
for i = 1:length(cell_idxs)
    spk_times = spikes{cell_idxs(i)} - trigger;
    spk_idxs = find(spk_times >= start & spk_times <= stop);
    spk_times = spk_times(spk_idxs);
    for j = 1:length(spk_times)
        bin_spikes(i,:,k) = bin_spikes(i,:,k) + (spk_times(j) >= bin_spikes(1:end-1) && spk_times(j) < bin_spikes(2:end));
    end
end
end

p_s = mean(mean(bin_spikes,3),2) ./ mean(bin_spikes(:));
I = mean(p_s .* log2(p_s + eps));

