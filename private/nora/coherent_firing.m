clear 
datarun_class = load_data('/Volumes/Analysis/2015-05-27-11/data001-data002/data001/data001');
datarun_class = load_params(datarun_class);
datapath = '/Volumes/Analysis/2015-05-27-11/data001-data002/data002/data002';
datarun= load_data(datapath);
datarun = load_neurons(datarun);

%% Find block starts
triggers = datarun.triggers;
trigger_diff = diff(triggers);
trigger_diff = abs(trigger_diff - median(trigger_diff));
repeat_starts = triggers([true; trigger_diff > 0.1]);
block_starts = [0; triggers([false; trigger_diff > 2])];

%% Find the different stim blocks
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

% load the test movie
disp('making test movie')
total_frames = 20*120;
testmovie = zeros(80, 40, total_frames, 'uint8');
total_images = total_frames/120;
for i = 1:total_images
    k = i; % using movie A for testing
    load(['/Volumes/Data/Stimuli/movies/eye-movement/current_movies/NSbrownian/matfiles/movie_chunk_' num2str(k) '.mat'])
    bin_movie = imresize(movie, 0.25);
    idx = (1:120) + 120*(i-1);
    testmovie(:,:,idx) = uint8(bin_movie);
    if ~mod(i, 50)
        disp(i)
    end
end
%%
figure; hold on
imagesc(testmovie(:,:,1081)'); axis image; colormap gray
cells = get_cell_indices(datarun_class, 'On Parasol');
for cell = cells
    center = datarun_class.vision.sta_fits{cell}.mean;
    text(center(1), 40-center(2), num2str(cell), 'Color', 'r');
end

%%
pair = [382 368];
%pair = [90 130];
for i_pair = 1:2
    tspikes{i_pair} = [];
    spikes = datarun.spikes{pair(i_pair)};
    testblocks = NSEM(1:2:end);
    testmovie_frames_per_block = 20*120;
    testmovie_seconds_per_block = testmovie_frames_per_block/120;
    for i = 1:length(testblocks)
        tspikes{i_pair} = [tspikes{i_pair}; spikes(spikes > testblocks(i)+9 & spikes < testblocks(i)+10)];
    end
end

[ccf, time] = compute_ccf(tspikes{1}, tspikes{2});
figure; plot(time, ccf)
