%% Load data
datapath = '2015-05-27-11/data001-data005-norefit/data002-from-data001_data002_data003_data004_data005/data002-from-data001_data002_data003_data004_data005';
datarun= load_data(datapath);
datarun = load_neurons(datarun);

% datapath = '2015-05-27-11/data001-data005-norefit/data001-from-data001_data002_data003_data004_data005/data001-from-data001_data002_data003_data004_data005';
% datarun_class= load_data(datapath);
% datarun_class = load_neurons(datarun_class);
% datarun_class = load_params(datarun_class);


%% Find block starts
triggers = datarun.triggers;
trigger_diff = diff(triggers);
trigger_diff = abs(trigger_diff - median(trigger_diff));
% hist(trigger_diff,1000)
% there appear to be a bunch of triggers around 0.2 and a bunch around 9
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

%%
% % get the right NSEM movie for the rasters
% disp('making test movie')
% testmovie_frames = 2400;
% testmovie_images = testmovie_frames/120;
% testmovie = zeros(320, 160, testmovie_frames, 'uint8');
% testmovie = uint8(testmovie);
% for i = 1:testmovie_images
%     load(['/Volumes/Data/Stimuli/movies/eye-movement/current_movies/NSbrownian_6000/matfiles/movie_chunk_' num2str(i) '.mat'])
%     idx = (1:120) + 120*(i-1);
%     testmovie(:,:,idx) = uint8(movie);
% end

%%
% load the fitting movie
disp('making fit movie')
total_frames = 30*7200;
fitmovie = zeros(80, 40, total_frames, 'uint8');
total_images = total_frames/120;
for i = 1:total_images
    k = i+3020; % using movie B
    load(['/Volumes/Data/Stimuli/movies/eye-movement/current_movies/NSbrownian_6000/matfiles/movie_chunk_' num2str(k) '.mat'])
    bin_movie = imresize(movie, 0.25);
    idx = (1:120) + 120*(i-1);
    fitmovie(:,:,idx) = uint8(bin_movie);
    if ~mod(i, 50)
        disp(i)
    end
end
disp('saving movie')
tic
save('/Volumes/Lab/Users/Nora/NSEM_Home/Stimuli/fitmovie_2015_05_27_11_data002.mat', 'fitmovie', '-v7.3')
toc

%% 
% load the test movie
disp('making test movie')
total_frames = 20*120;
testmovie = zeros(80, 40, total_frames, 'uint8');
total_images = total_frames/120;
for i = 1:total_images
    k = i; % using movie A for testing
    load(['/Volumes/Data/Stimuli/movies/eye-movement/current_movies/NSbrownian_6000/matfiles/movie_chunk_' num2str(k) '.mat'])
    bin_movie = imresize(movie, 0.25);
    idx = (1:120) + 120*(i-1);
    testmovie(:,:,idx) = uint8(bin_movie);
    if ~mod(i, 50)
        disp(i)
    end
end
disp('saving movie')
tic
save('/Volumes/Lab/Users/Nora/NSEM_Home/Stimuli/testmovie_2015_05_27_11_data002.mat', 'testmovie', '-v7.3')
toc
