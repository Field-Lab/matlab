clear

%% Load data
datapath = '2015-05-27-11/data001-data005-norefit/data002-from-data001_data002_data003_data004_data005/data002-from-data001_data002_data003_data004_data005';
datarun= load_data(datapath);
datarun = load_neurons(datarun);

datapath = '2015-05-27-11/data001-data005-norefit/data001-from-data001_data002_data003_data004_data005/data001-from-data001_data002_data003_data004_data005';
datarun_class= load_data(datapath);
datarun_class = load_neurons(datarun_class);
datarun_class = load_params(datarun_class);

%%
% disp('loading fit movie')
% load('/Volumes/Lab/Users/Nora/NSEM_Home/Stimuli/fitmovie_2015_05_27_11_data002.mat')
disp('loading test movie')
load('/Volumes/Lab/Users/Nora/NSEM_Home/Stimuli/testmovie_2015_05_27_11_data002.mat')
testmovie = permute(testmovie, [2 1 3]);
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
% Organize test spikes
testblocks = NSEM(1:2:end);
testmovie_frames_per_block = 20*120;
testmovie_seconds_per_block = testmovie_frames_per_block/120;

cell_type = {'On Parasol', 'Off Amacrine', 'Off Parasol'};
for i_cell_type = 1:3

cells = get_cell_indices(datarun_class, cell_type{i_cell_type});
n_cells = length(cells);

res{i_cell_type}.spikes = zeros(n_cells, testmovie_frames_per_block);
res{i_cell_type}.centers = zeros(n_cells, 2);

for i_cell = 1:n_cells
    
    disp(i_cell)
    spikes = datarun.spikes{cells(i_cell)};
    % Concatenate fit spikes
%     fitblocks = NSEM(2:2:end);
%     fitmovie_frames_per_block = 7200;
%     fitmovie_seconds_per_block = fitmovie_frames_per_block/120;
%     concat_spikes = [];
%     start = 0;
%     for i = 1:length(fitblocks)
%         block_spikes = spikes(spikes > fitblocks(i) & spikes < fitblocks(i)+fitmovie_seconds_per_block);
%         concat_spikes = [concat_spikes; block_spikes-fitblocks(i)+start];
%         start = start + fitmovie_seconds_per_block;
%     end

    spikes_concat{i_cell_type}{i_cell} = [];
    for i = 1:length(testblocks)
        testspikes{i_cell, i} = spikes(spikes > testblocks(i) & spikes < testblocks(i)+testmovie_seconds_per_block) - testblocks(i);
        spikes_concat{i_cell_type}{i_cell} = [spikes_concat{i_cell_type}{i_cell}; testspikes{i_cell, i}];
    end
%      spikes_frame = floor(spikes_concat{i_cell} * 120);
%      for i_frame = 1:testmovie_frames_per_block
%         res{i_cell_type}.spikes(i_cell, i_frame) = sum(spikes_frame == i_frame); 
%      end

    res{i_cell_type}.centers(i_cell,:) = datarun_class.vision.sta_fits{cells(i_cell)}.mean;
    %center(2) = 40 - center(2);

end
end

%%
tic; res_spikes_plot(testmovie, res, ['/Users/Nora/Desktop/' cell_type{i_cell_type} '_movie.avi']); toc

%%
default_colors = get(gca,'ColorOrder');
for i_cell_type = 1:3
    subplot(3,1,i_cell_type)
    for i_cell = 1:size(res{i_cell_type}.centers,1)
        hold on; plot(spikes_concat{i_cell_type}{i_cell},res{i_cell_type}.centers(i_cell,1)*ones(size(spikes_concat{i_cell_type}{i_cell})), '.k');
    end
    title(cell_type{i_cell_type})
end
set(gcf, 'Position', [100 100 800 650])
for i_saccade = 0:19
    for i_cell_type = 1:3
    subplot(3,1,i_cell_type)
    xlim([i_saccade i_saccade+0.25])
    end
    %suptitle(['Saccade ' num2str(i_saccade)])
    %hold on; subplot(2,1,2); plot(spikes, centers(2)*ones(length(spikes), 1), '.k')
    %xlim([i_saccade i_saccade+0.25])
    exportfig(gcf, ['/Users/Nora/Desktop/waves/Saccade' num2str(i_saccade)], 'Bounds', 'loose', 'Color', 'rgb', 'Renderer', 'opengl');
end
%%
