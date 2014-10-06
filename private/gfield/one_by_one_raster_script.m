
temp_path = '/snle/lab/Experiments/Array/Analysis/2008-05-13-3/data008-mg/data008/data008';
datarun = load_data(temp_path);
datarun = load_sta(datarun, 'load_sta', []);
datarun = load_neurons(datarun);


cycle_duration = 90; %seconds


%%% RASTER
cell_ID = 2404; %off midget
cell_ID = 1684 % off parasol
cell_index = get_cell_indices(datarun, cell_ID);
triggers_per_cycle = 109;
begin_time = 0;
end_time = 1800;
trigger_begin = ((begin_time ./ cycle_duration) * triggers_per_cycle) +1;
trigger_end = (end_time ./ cycle_duration) .* triggers_per_cycle;

cycle_trigger_indices = trigger_begin:triggers_per_cycle:trigger_end;
num_cycles = length(cycle_trigger_indices)-1;

cycle_triggers = datarun.triggers(cycle_trigger_indices);
temp_spikes = datarun.spikes{cell_index};

[raster_times, raster_marks] = spike_raster(temp_spikes, cycle_triggers, datarun.triggers(1), 'plot_raster', true, 'raster_size', 16);

