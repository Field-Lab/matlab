
clear
reg{1} = 'data003';
reg{2} = 'data007';
reg{3} = 'data011';

classification = 'data001';

Analysis_Path = '/Volumes/Analysis/2016-01-05-0/map-from-data001/';
fig_save = '/Users/Nora/Desktop/Fig_Output';
mkdir(fig_save);

set(0, 'defaultFigurePaperPositionMode', 'auto')
set(0, 'defaultFigurePaperOrientation', 'landscape')
set(0, 'defaultFigureUnits', 'inches')
set(0, 'defaultFigurePosition', [2 2 8 6])

% classification run
class_datarun = load_data(['/Volumes/Analysis/2016-01-05-0/streamed/' classification '/' classification]);
class_datarun = load_params(class_datarun);
class_datarun = load_neurons(class_datarun);

%% load NSinterval: only for GS!

i_chunk = 1;
load(['/Users/Nora/Desktop/Data/matfiles/movie_chunk_' num2str(i_chunk) '.mat']);
NSmovie = movie_chunk;
i_chunk = 2;

% mask movie
while size(NSmovie,3) < 1200
    load(['/Users/Nora/Desktop/Data/matfiles/movie_chunk_' num2str(i_chunk) '.mat']);
    NSmovie = cat(3,NSmovie, movie_chunk);
    if size(NSmovie,3) > 1200
        NSmovie = NSmovie(:,:,1:1200);
    end
    i_chunk = i_chunk + 1;
end
NSmovie = permute(NSmovie, [2 1 3]);


%%
i = 3;
datarun = load_data([ Analysis_Path reg{i} '/' reg{i}]);
datarun = load_neurons(datarun);
data = interleaved_data_prep(datarun, 1100, 30, 'cell_spec', 'Off Amacrine', 'datarun_class', class_datarun, 'visual_check', 0);
cell_indices = get_cell_indices(class_datarun, 'Off Amacrine');
n_cells = length(cell_indices);
res.spikes = zeros(n_cells, 1200);
for i_cell = 1:n_cells
    res.spikes(i_cell, :) = IDP_plot_PSTH(data, i_cell);
    res.centers(i_cell,:) = 8*class_datarun.vision.sta_fits{cell_indices(i_cell)}.mean;
end
res_spikes_plot(NSmovie, res, '/Users/Nora/Desktop/amacrine_movie3.avi')

%%

% get image transition times
% mask movie
i_chunk = 1;
image_transitions = [1];
while image_transitions(end) < 1200
    load(['/Users/Nora/Desktop/Data/matfiles/movie_chunk_' num2str(i_chunk) '.mat']);
    image_transitions = [image_transitions image_transitions(end)+size(movie_chunk,3)];
    i_chunk = i_chunk + 1;
end
image_transitions = image_transitions(1:(end));

%%
figure; hold on
slope = zeros(20, 2);

for i_cell = 1:length(cell_indices)
    spikes = cell2mat(data.testspikes(:,i_cell));
    centers = flip(class_datarun.vision.sta_fits{cell_indices(i_cell)}.mean);
    hold on; subplot(2,1,1); plot(spikes, centers(1)*ones(length(spikes), 1), '.k')
    hold on; subplot(2,1,2); plot(spikes, centers(2)*ones(length(spikes), 1), '.k')
end
for i_saccade = image_transitions/120;
    subplot(2,1,1);
    title(['Saccade ' num2str(i_saccade)])
    xlim([i_saccade i_saccade+0.25])
    subplot(2,1,2);
    xlim([i_saccade i_saccade+0.25]) 
    pause()
end


