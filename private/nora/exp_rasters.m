datapath = '2015-04-14-1/data000-data004-norefit/data003-from-data000_data001_data002_data003_data004/data003-from-data000_data001_data002_data003_data004';
datarun = load_data(datapath);
datarun = load_neurons(datarun);

%%
triggers = datarun.triggers;
trigger_diff = diff(triggers);
trigger_diff = abs(trigger_diff - median(trigger_diff));

%%
hist(trigger_diff,1000)
% there appear to be a bunch of triggers around 0.2 and a bunch around 9

%%
repeat_starts = triggers([true; trigger_diff > 0.1]);
block_starts = triggers([true; trigger_diff > 2]);


%%
all_rasters = repeat_starts(1:2:end);
rasters{1} = all_rasters([1:4,11:14, 21:24, 31:34, 41:44]);
rasters{2} = all_rasters(5+[1:4,11:14, 21:24, 31:34, 41:44]);

%%
raster_length = 30;
cell = 72;
stim = 2;
hold on
for i = 1:length(rasters{stim})
    trial_spikes{i} = datarun.spikes{cell}(datarun.spikes{cell} > rasters{stim}(i) & datarun.spikes{cell} < (rasters{stim}(i) + raster_length))-rasters{stim}(i);
    plot(trial_spikes{i},i*ones(length(trial_spikes{i}),1),'.k')
end

%%

hold on
plot(repeat_starts-1, 0.3*ones(length(repeat_starts),1),'*')

%%
plot(triggers, [0; trigger_diff])
hold on
plot(rasters{2}, 0.1*ones(length(rasters{2})),'*')
