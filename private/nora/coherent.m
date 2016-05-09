%% Dataruns: run first for everything

clear
reg{1} = 'data003'; 
reg{2} = 'data007'; 
reg{3} = 'data011';

classification = 'data001';

Analysis_Path = '/Volumes/Analysis/2016-01-05-0/map-from-data001/';
fig_save = '/Users/Nora/Desktop/Fig_Output';
%mkdir(fig_save);

set(0, 'defaultFigurePaperPositionMode', 'auto')
set(0, 'defaultFigurePaperOrientation', 'landscape')
set(0, 'defaultFigureUnits', 'inches')
set(0, 'defaultFigurePosition', [2 2 15 4])

% classification run
class_datarun = load_data(['/Volumes/Analysis/2016-01-05-0/streamed/' classification '/' classification]);
class_datarun = load_params(class_datarun);
class_datarun = load_neurons(class_datarun);

%% check stability
maskin = 1;
datarun = load_data([ Analysis_Path reg{maskin} '/' reg{maskin}]);
datarun = load_neurons(datarun);

%% load NSinterval: only for GS!
i_chunk = 1;
load(['/Volumes/Data/Stimuli/movies/eye-movement/current_movies/NSinterval/matfiles/movie_chunk_' num2str(i_chunk) '.mat']);
NSmovie = movie_chunk;
interval(i_chunk) = size(movie_chunk, 3);

i_chunk = 2;

% mask movie
while size(NSmovie,3) < 1200
    load(['/Volumes/Data/Stimuli/movies/eye-movement/current_movies/NSinterval/matfiles/movie_chunk_' num2str(i_chunk) '.mat']);
    interval(i_chunk) = size(movie_chunk, 3);
    
    NSmovie = cat(3,NSmovie, movie_chunk);
    if size(NSmovie,3) > 1200
        NSmovie = NSmovie(:,:,1:1200);
    end
    i_chunk = i_chunk + 1;
    %% SECTION TITLE
    % DESCRIPTIVE TEXT
end
clear movie_chunk
interval_frame = cumsum(interval);
interval_time = interval_frame * (1/119.5);

%%
fix = 25;
imagesc(NSmovie(:,:,interval_frame(fix))'); axis image; colormap gray
cell_ids = get_cell_ids(class_datarun, 'On Parasol');
cell_idx = get_cell_indices(class_datarun, 'On Parasol');
n_cells = length(cell_ids);
reg_data = interleaved_data_prep(datarun, 1100, 30, 'cell_spec', 'On Parasol', 'datarun_class', class_datarun, 'visual_check', 0);

hold on; 
for i_cell = 1:n_cells
    center = 8*class_datarun.vision.sta_fits{cell_idx(i_cell)}.mean;
    text(center(1), 160-center(2), num2str(i_cell), 'FontSize', 14, 'Color', 'r')
end

%%
test_pair = [37 19];
trial_length = 1;
spikes{1} = [];
spikes{2} = [];
for i_trial = 1:size(reg_data.testspikes, 1)
    for pair = 1:2
        trial_spikes = reg_data.testspikes{i_trial, test_pair(pair)};
        trial_spikes_fix = trial_spikes((trial_spikes > interval_time(fix-1)) & (trial_spikes < interval_time(fix)))+trial_length*i_trial;
        spikes{pair} = [spikes{pair}; trial_spikes_fix];
    end
end

figure(2); 
[ccf, time] = compute_ccf(spikes{1}, spikes{2});
plot(time, ccf)

