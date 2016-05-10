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

WN4 = [];
WN8 = [];
NSEM = [];

for i = 1:length(block_starts)-1
    repeats_within_block = repeat_starts(repeat_starts > block_starts(i) & repeat_starts < block_starts(i+1));
    repeats_within_block = repeats_within_block(1:20);
    if mod(i,3) == 1
        WN4 = [WN4; repeats_within_block];
    elseif mod(i,3) == 2
        WN8 = [WN8; repeats_within_block];
    else
        NSEM = [NSEM; repeats_within_block];
    end
end

%%
amacrines = get_cell_indices(datarun_class, 'Off Amacrine');

for i =8 % 1:length(amacrines)
    % Organize test spikes
%     testblocks = NSEM(1:2:end);
%     for i = 1:length(testblocks)
%         tspikes{i} = spikes(spikes > testblocks(i) & spikes < testblocks(i)+testmovie_seconds_per_block) - testblocks(i);
%     end
%     
    cell = amacrines(i);
    starts = NSEM(1:2:end);
    starts = [starts(1:10); starts(13:20); starts(23:end)];
    spikes = get_raster(datarun.spikes{cell}, starts, 'stop', 10, 'start', 2);
    %     a=text(10.1, 5, '4 stixel WN', 'Color', 'r', 'FontSize', 8);
    %     b=text(10.1, 15, '8 stixel WN', 'Color', 'r', 'FontSize', 8);
    %     c=text(10.1, 25, 'NSEM', 'Color', 'r', 'FontSize', 8);
    % pause()
    set(gcf, 'Position', [100 100 800 250])
    %xlabel('Time (seconds)')
    %ylabel('Trials')
    axis off
    
end
