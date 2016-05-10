clear
Analysis_Path = '/Volumes/Analysis/2016-02-17-1/data022-data028';
datarun_class = load_data([Analysis_Path '/data024/data024'], struct('load_neurons', 0, 'load_params', 1));

%% only change things here
cell_type = 'Off Parasol';
cells = get_cell_ids(datarun_class, cell_type); % cell ids to fit
cids = get_cell_indices(datarun_class, cell_type); % cell ids to fit

%% BW
% test_data = 'data027';
% test_datarun = load_data([Analysis_Path '/' test_data '/' test_data], struct('load_neurons', 1, 'load_params', 1));
% repeats = interleaved_data_prep(test_datarun, 1100, 30, 'cell_spec', cells,'visual_check', 0, 'stimulus_name', 'BW-8-1', 'seed', 22222);

%% NSEM
test_data = 'data022';
test_datarun = load_data([Analysis_Path '/' test_data '/' test_data], struct('load_neurons', 1, 'load_params', 1));
repeats = interleaved_data_prep(test_datarun, 1100, 30,'cell_spec', cells,'visual_check', 0);

set(gcf, 'Position', [100 100 800 250])
for i_cell = 1:size(repeats.testspikes,2)
        spikes_concat = cell2mat(repeats.testspikes(:,i_cell));
        centers = flip(datarun_class.vision.sta_fits{cids(i_cell)}.mean);
        hold on; plot(spikes_concat,centers(1)*ones(size(spikes_concat)), '.k');
end
ylim([0 35])
for i_saccade = 0:9
    xlim([i_saccade i_saccade+0.25])
    pause()
    %exportfig(gcf, ['/Users/Nora/Desktop/waves/Saccade' num2str(i_saccade)], 'Bounds', 'loose', 'Color', 'rgb', 'Renderer', 'opengl');
end