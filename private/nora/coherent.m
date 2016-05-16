%% Dataruns: run first for everything

%{
clear
reg{1} = 'data005'; 
reg{2} = 'data007'; 
reg{3} = 'data011';
classification = 'data002';
Analysis_Path = '/Volumes/Analysis/2015-12-18-2/data002-data015/';
class = [Analysis_Path classification '/' classification];
%}

piece = '2016-02-17-1';
Analysis_Path = ['/Volumes/Analysis/' piece '/mVision/'];
class = ['/Volumes/Analysis/' piece '/streamed/data019/data019'];
reg{1} = 'data022';
long = 'data025';
%}

set(0, 'defaultFigurePaperPositionMode', 'auto')
set(0, 'defaultFigurePaperOrientation', 'landscape')
set(0, 'defaultFigureUnits', 'inches')
set(0, 'defaultFigurePosition', [2 2 15 4])

% classification run
class_datarun = load_data(class);
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
end
clear movie_chunk
interval_frame = cumsum(interval);
interval_time = interval_frame * (1/119.5);

%%
cell_ids = get_cell_ids(class_datarun, 'On Parasol');
%cell_ids = [cell_ids get_cell_ids(class_datarun, 'Off Parasol')];
cell_idx = get_cell_indices(class_datarun, 'On Parasol');
ref_data = '/Volumes/Analysis/2016-02-17-1/streamed/data019';
new_data = '/Volumes/Analysis/2016-02-17-1/mVision/data022';
EI_map = crossIdentifyNeuronIDs_NB(ref_data, new_data, cell_ids, [], 0);
cell_ids = EI_map(:,2);
cell_ids = cell_ids(~isnan(cell_ids));
n_cells = length(cell_ids);
reg_data = interleaved_data_prep(datarun, 1100, 30, 'cell_spec', cell_ids, 'visual_check', 0);
%%
test_pair = [1 3]+18; %10, 12
center(1,:) = 4*class_datarun.vision.sta_fits{cell_idx(test_pair(1))}.mean;
center(2,:) = 4*class_datarun.vision.sta_fits{cell_idx(test_pair(2))}.mean;
spikes{1} = [];
spikes{2} = [];
trial_length = 10; 

%
for fix = 2:25
    bin= im2bw(NSmovie(:,:,interval_frame(fix))'/255, 0.25);
    temp = improfile(bin, [center(1,1), center(2,1)],[160-center(1,2) 160-center(2,2)] );
   %imagesc(bin);
   %hold on; plot([center(1,1), center(2,1)],[160-center(1,2) 160-center(2,2)]); pause()
    if std(temp) == 0
        disp('sameobj')
        for i_trial = 1:size(reg_data.testspikes, 1)
            for pair = 1:2
                trial_spikes = reg_data.testspikes{i_trial, test_pair(pair)};
                trial_spikes_fix = trial_spikes((trial_spikes > interval_time(fix-1)) & (trial_spikes < interval_time(fix)))+trial_length*i_trial;
                spikes{pair} = [spikes{pair}; trial_spikes_fix];
            end
        end
    end
end

%
figure(2); 
[ccf, time] = compute_ccf(sort(spikes{1}), sort(spikes{2}));
plot(time, ccf)

%
