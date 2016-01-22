% %% Load data
datapath = '2015-05-27-11/data001-data005-norefit/data002-from-data001_data002_data003_data004_data005/data002-from-data001_data002_data003_data004_data005';
datarun= load_data(datapath);
datarun = load_neurons(datarun);

datapath = '2015-05-27-11/data001-data005-norefit/data001-from-data001_data002_data003_data004_data005/data001-from-data001_data002_data003_data004_data005';
datarun_class= load_data(datapath);
datarun_class = load_neurons(datarun_class);
datarun_class = load_params(datarun_class);


%% Find block starts
% Finding block starts
triggers = datarun.triggers;
trigger_diff = diff(triggers);
trigger_diff = abs(trigger_diff - median(trigger_diff));
% hist(trigger_diff,1000)
% there appear to be a bunch of triggers around 0.2 and a bunch around 9
repeat_starts = triggers([true; trigger_diff > 0.1]);
block_starts = [0; triggers([false; trigger_diff > 2])];

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

clear WN8 WN4 repeats_within_block repeat_starts block_starts

%% Load up cell info
cells1 = get_cell_indices(datarun_class, 'On Parasol');
cells2 = get_cell_indices(datarun_class, 'Off Amacrine');
hold on
for i_cell1 = 1:length(cells1)
    disp(i_cell1)
    for i_cell2 = 1:length(cells2)
        [ccf{i_cell1, i_cell2}, time{i_cell1, i_cell2}] = correlation_shift(datarun.spikes{cells1(i_cell1)},datarun.spikes{cells2(i_cell2)});
        % plot(time{i_cell1, i_cell2}, ccf{i_cell1, i_cell2})
        % pause(0.01)
    end
end
save('CCF.mat', 'ccf')

%%
%distance2 = zeros(7,25);
% CCFpeak = zeros(53*29, 3);
%count = 0;
hold on
for i_cell1 = 1:53
    for i_cell2 = 1:29
        %distance2(i_cell1, i_cell2) = 
        %count = count + 1;
        % CCFpeak(count, 1) = norm(datarun_class.vision.sta_fits{cells1(i_cell1)}.mean - datarun_class.vision.sta_fits{cells2(i_cell2)}.mean);
        % [CCFpeak(count, 2), CCFpeak(count, 3)] = max(abs(ccf{i_cell1, i_cell2}));
        plot(time{i_cell1, i_cell2}, ccf{i_cell1, i_cell2})
        % figure(1); hold on;  plot(CCFpeak(count, 1), CCFpeak(count, 2), '*', 'Color', default_colors(i_cell1, :))
        % figure(2); hold on; plot(CCFpeak(count, 1), CCFpeak(count, 3), '*', 'Color', default_colors(i_cell1, :))
        % figure(3); hold on; plot(CCFpeak(count, 2), CCFpeak(count, 3), '*', 'Color', default_colors(i_cell1, :))
    end
end