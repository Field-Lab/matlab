% Script for generating absolute threshold plots for rod input to SBC paper
set(0,'DefaultAxesFontSize', 16, 'DefaultAxesFontName', 'Helvetica')
set(0, 'DefaultAxesColorOrder', [1 0 0; 0 1 0; 0 0 1; 0 0 0]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data loading and mapping

% set data paths and array type
% note, these data are suspected to be misleading because there may have
% been a light leak from the monitor into the faraday cage around the rig.
master_data_path = '/snle/lab/Experiments/Array/Analysis/2008-11-12-2/data000/data000';
slave_data_path = '/snle/lab/Experiments/Array/Analysis/2008-11-12-2/data006/data006-manual/data006-manual';
array_type = 519;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
master_data_path = '/snle/lab/Experiments/Array/Analysis/2008-12-12-1/data006/data006';
slave_data_path = '/snle/lab/Experiments/Array/Analysis/2008-12-12-1/data014/data014';  %NDF 5.0
slave_data_path = '/snle/lab/Experiments/Array/Analysis/2008-12-12-1/data015/data015';  %NDF 5.6
slave_data_path = '/snle/lab/Experiments/Array/Analysis/2008-12-12-1/data012/data012';  %NDF 6.0
slave_data_path = '/snle/lab/Experiments/Array/Analysis/2008-12-12-1/data017/data017';  %NDF 6.0
slave_data_path = '/snle/lab/Experiments/Array/Analysis/2008-12-12-1/data017/data017-manual/data017-manual';  %NDF 6.0
slave_data_path = '/Analysis/gfield/2008-12-12-1/data012-data018/data012-data018'; %NDF inf

slave_data_path = '/snle/lab/Experiments/Array/Analysis/2008-12-12-1/data016/data016';  %NDF 6.3
slave_data_path = '/snle/lab/Experiments/Array/Analysis/2008-12-12-1/data016/data016';  %NDF 6.6
slave_data_path = '/snle/lab/Experiments/Array/Analysis/2008-12-12-1/data013/data013';  %NDF 6.9
slave_data_path = '/snle/lab/Experiments/Array/Analysis/2008-12-12-1/data020/data020';   %NDF 6.3

array_type = 512;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load master data
clear temp_datarun
temp_datarun = load_data(master_data_path);
temp_datarun = load_params(temp_datarun,struct('verbose',1));  
temp_datarun = load_ei(temp_datarun, 'all', struct('array_type', array_type));
datarun{1} = temp_datarun;

% load slave data
clear temp_datarun
temp_datarun = load_data(slave_data_path);
temp_datarun = load_neurons(temp_datarun);
temp_datarun = load_ei(temp_datarun, 'all',struct('array_type', array_type));
datarun{2} = temp_datarun;

% EI mapping from RF run to control condition for SBCs
threshold = 0.9435; % correlaton threhold for mapping
cell_type = {5};
master_cell_type = cell_type;
slave_cell_type = 'all';
[datarunA, datarunB] = map_gdf(datarun{1}, datarun{2}, 'corr_threshold', threshold, 'master_cell_type', master_cell_type, 'slave_cell_type', slave_cell_type, 'verbose', true, 'troubleshoot', true);

% EI mapping for ON parasols
threshold = 0.95; % correlaton threhold for mapping
cell_type = {1};
master_cell_type = cell_type;
slave_cell_type = 'all';
[datarunA, datarunB] = map_gdf(datarunA, datarunB, 'corr_threshold', threshold, 'master_cell_type', master_cell_type, 'slave_cell_type', slave_cell_type, 'verbose', true, 'troubleshoot', true);

% EI mapping for OFF parasols
threshold = 0.95; % correlaton threhold for mapping
cell_type = {2};
master_cell_type = cell_type;
slave_cell_type = 'all';
[datarunA, datarunB] = map_gdf(datarunA, datarunB, 'corr_threshold', threshold, 'master_cell_type', master_cell_type, 'slave_cell_type', slave_cell_type, 'verbose', true, 'troubleshoot', true);

% EI mapping for ON midgets
threshold = 0.95; % correlaton threhold for mapping
cell_type = {3};
master_cell_type = cell_type;
slave_cell_type = 'all';
[datarunA, datarunB] = map_gdf(datarunA, datarunB, 'corr_threshold', threshold, 'master_cell_type', master_cell_type, 'slave_cell_type', slave_cell_type, 'verbose', true, 'troubleshoot', true);

% EI mapping for OFF midgets
threshold = 0.95; % correlaton threhold for mapping
cell_type = {4};
master_cell_type = cell_type;
slave_cell_type = 'all';
[datarunA, datarunB] = map_gdf(datarunA, datarunB, 'corr_threshold', threshold, 'master_cell_type', master_cell_type, 'slave_cell_type', slave_cell_type, 'verbose', true, 'troubleshoot', true);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% extract raster
% SBC 5194
% on-parasol 5228
temp_datarun = datarunB;
cell_number = 18;
cell_type = {3};
cell_ID = temp_datarun.cell_types{cell_type{1}}.cell_ids(cell_number)
cell_index = get_cell_indices(temp_datarun, cell_ID);

cycle_duration = 20; %seconds
bin_interval = 0.2; %seconds
triggers_per_cycle = 4;
begin_time = 0;
%end_time = 400;
end_time = temp_datarun.duration;
trigger_begin = ((begin_time ./ cycle_duration) * triggers_per_cycle) +1;
trigger_end = (end_time ./ cycle_duration) * triggers_per_cycle;

hist_bins = 0:bin_interval:cycle_duration;
cycle_trigger_indices = trigger_begin:triggers_per_cycle:trigger_end;
%cycle_trigger_indices = 1:triggers_per_cycle:length(temp_datarun.triggers);
num_cycles = length(cycle_trigger_indices)-1;
num_bins = length(hist_bins);

cycle_triggers = temp_datarun.triggers(cycle_trigger_indices);
spikes = temp_datarun.spikes{cell_index};

[raster_times, raster_marks] = spike_raster(spikes, cycle_triggers, temp_datarun.triggers(1), 'plot_raster', true, 'raster_size', 16);
print(1, '/snle/home/gfield/Desktop/off-par-raster','-dpdf')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate PSTHs for all cells

control_cell_indices = get_cell_indices(datarunB, cell_type);
cell_number = length(control_cell_indices);

for cll = 1:cell_number;
    temp_control_spikes = datarunB.spikes{control_cell_indices(cll)};
    temp_datarun = datarunB;
    hist_matrix = zeros(num_cycles, num_bins);
    for cntr = 1:num_cycles
        temp_index = cycle_trigger_indices(cntr);
        temp_index_2 = cycle_trigger_indices(cntr+1);
        start_time_interval = temp_datarun.triggers(temp_index) - temp_datarun.triggers(1);
        end_time_interval = temp_datarun.triggers(temp_index_2) - temp_datarun.triggers(1);

        above_indices = find(temp_control_spikes > start_time_interval);
        below_indices = find(temp_control_spikes < end_time_interval);
        temp_spike_indices = intersect(above_indices, below_indices);

        temp_cycle_spike_times = temp_control_spikes(temp_spike_indices);

        temp_hist_bins = hist_bins + start_time_interval;
        temp_hist = hist(temp_cycle_spike_times, temp_hist_bins) ./ bin_interval;

        hist_matrix(cntr,:) = temp_hist;
    end
    temp_cycle_hist = sum(hist_matrix, 1) ./ num_cycles;
    control_cycle_hist(cll,:) = temp_cycle_hist;
end

% Generate stimulus trace
trigger_one = datarunB.triggers(1);
trigger_two = datarunB.triggers(2);
trigger_three = datarunB.triggers(3);
trigger_four = datarunB.triggers(4);
stim_trace_matrix = [0 trigger_one-0.001 trigger_one trigger_two-0.001 trigger_two trigger_three-0.001 trigger_three trigger_four-0.001 trigger_four cycle_duration;...
                     -2 -2 -1 -1 -2 -2 -3 -3 -2 -2];
 
%plot histogram across cells
cell_ids = datarunB.cell_types{cell_type{1}}.cell_ids;
%for cll = 1:length(cell_ids)
for cll = 18
%for cll = find(cell_ids == cell_of_interest)
    figure(cll)
    clf
    bar(hist_bins, control_cycle_hist(cll, :), 'k')
    ylabel('Hz')
    hold on
    % plots stimulus trace
    plot(stim_trace_matrix(1,:), stim_trace_matrix(2,:), 'k') 
    axis([0 20 -4 60])
    hold off
end
print(cll, '/snle/home/gfield/Desktop/off-par-hist','-dpdf')


on_midget_example_list = [48 44 42 31];
on_midget_exclusion_list = [46 38 32 25 19 12 5 3 2];

off_midget_example_list = [93 119 116 112];
off_midget_exclusion_list = [125 118 109 107 106 105 99 91 52 51 50 47 46 44 42 41,...
                            40 35 30 28 27 26 25 20 12 1];
                        
off_par_example_list = [20 18 29 23 22];                        
off_par_exclusion_list = [33 31 27 25 14 12 9 7 6 3 2];                      

% exclude bad cells from mean and std
temp_cell_list = 1:length(cell_ids);
if cell_type{1} == 4
    exclude_list = off_midget_exclusion_list;
elseif cell_type{1} == 3
    exclude_list = on_midget_exclusion_list;
elseif cell_type{1} == 2
    exclude_list = off_par_exclusion_list;
end
cell_list = setdiff(temp_cell_list, exclude_list);
                        
% average histograms across cells
mean_control_cycle_hist = mean(control_cycle_hist(cell_list,:), 1);
mean_control_cycle_hist(1) = mean_control_cycle_hist(3);
mean_control_cycle_hist(101) = mean_control_cycle_hist(98);

figure(200)
clf
hold on
plot(hist_bins, mean_control_cycle_hist, 'k')
axis([0 cycle_duration 0 20])
hold off
print(200, '/snle/home/gfield/Desktop/on-midget-mean-hist','-dpdf')

temp_index = find(datarunB.cell_types{5}.cell_ids == cell_of_interest);
cell_of_interest_master = datarunA.cell_types{5}.cell_ids(temp_index);
cell_of_interest_master_index = get_cell_indices(datarunA, cell_of_interest_master);
datarunA = load_sta(datarunA, 'load_sta', {5});
datarunA = set_polarities(datarunA, 'cell_specs', {{5}}, 'polarities', [-1]);

% plot spatial STA
datarunA = get_sta_summaries(datarunA, {5});
plot_rf_portraits(datarunA, cell_of_interest_master, 'plot_radius', 8, 'scale_factor', 10);
print(1, '/snle/home/gfield/Desktop/SBC-rf','-dpdf')

% plot STA time course
plot(datarunA.stas.time_courses{cell_of_interest_master_index});
axis([10 30 -0.06 0.06])
print(1, '/snle/home/gfield/Desktop/SBC-tc','-dpdf')




