datapath = '2015-05-27-11/data001-data005-norefit/data002-from-data001_data002_data003_data004_data005/data002-from-data001_data002_data003_data004_data005';
datarun= load_data(datapath);
datarun = load_neurons(datarun);

datapath = '2015-05-27-11/data001-data005-norefit/data001-from-data001_data002_data003_data004_data005/data001-from-data001_data002_data003_data004_data005';
datarun_class= load_data(datapath);
datarun_class = load_neurons(datarun_class);
datarun_class = load_params(datarun_class);
%%
triggers = datarun.triggers;
trigger_diff = diff(triggers);
trigger_diff = abs(trigger_diff - median(trigger_diff));

%%
% hist(trigger_diff,1000)
% there appear to be a bunch of triggers around 0.2 and a bunch around 9

%%
repeat_starts = triggers([true; trigger_diff > 0.1]);
block_starts = triggers([true; trigger_diff > 2]);
all_rasters = repeat_starts(1:2:end);

%%
amacrines = get_cell_indices(datarun_class, 'Off Amacrine');

for i = 1:length(amacrines)
    cell = amacrines(i);
    spikes = get_raster(datarun.spikes{cell}, all_rasters, 'stop', 10);
    a=text(10.1, 5, '4 stixel WN', 'Color', 'r', 'FontSize', 8);
    b=text(10.1, 15, '8 stixel WN', 'Color', 'r', 'FontSize', 8);
    c=text(10.1, 25, 'NSEM', 'Color', 'r', 'FontSize', 8);
    pause()
end
