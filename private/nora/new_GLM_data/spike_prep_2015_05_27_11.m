%% Load data
datapath = '2015-05-27-11/data001-data005-norefit/data002-from-data001_data002_data003_data004_data005/data002-from-data001_data002_data003_data004_data005';
datarun= load_data(datapath);
datarun = load_neurons(datarun);

datapath = '2015-05-27-11/data001-data005-norefit/data001-from-data001_data002_data003_data004_data005/data001-from-data001_data002_data003_data004_data005';
datarun_class= load_data(datapath);
datarun_class = load_neurons(datarun_class);
datarun_class = load_params(datarun_class);

%%
disp('loading fit movie')
load('/Volumes/Lab/Users/Nora/NSEM_Home/Stimuli/fitmovie_2015_05_27_11_data002.mat')
disp('loading test movie')
load('/Volumes/Lab/Users/Nora/NSEM_Home/Stimuli/testmovie_2015_05_27_11_data002.mat')

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

% cells{1} = [10 14 45 4];
% cells{2} = [285 268 339 5386];
% cells{1} = [14 4];

% cells{1} = [206 207 183 148];

% On par
cells{1} = 41;

%% Load up cell info

for i_cell = 1%:length(cells)
    
    disp(i_cell)
    cell = cells{i_cell}(1);
    spikes = datarun.spikes{cell};
    % Concatenate fit spikes
    fitblocks = NSEM(2:2:end);
    fitmovie_frames_per_block = 7200;
    fitmovie_seconds_per_block = fitmovie_frames_per_block/120;
    concat_spikes = [];
    start = 0;
    for i = 1:length(fitblocks)
        block_spikes = spikes(spikes > fitblocks(i) & spikes < fitblocks(i)+fitmovie_seconds_per_block);
        concat_spikes = [concat_spikes; block_spikes-fitblocks(i)+start];
        start = start + fitmovie_seconds_per_block;
    end
    % Organize test spikes
    testblocks = NSEM(1:2:end);
    testmovie_frames_per_block = 20*120;
    testmovie_seconds_per_block = testmovie_frames_per_block/120;
    for i = 1:length(testblocks)
        tspikes{i} = spikes(spikes > testblocks(i) & spikes < testblocks(i)+testmovie_seconds_per_block) - testblocks(i);
    end
    
    center = datarun_class.vision.sta_fits{cell}.mean;
    center(2) = 40 - center(2);
    
    cell_pairs = cells{i_cell}(2:end);
    for cell = 1:length(cell_pairs)
        spikes = datarun.spikes{cell};
        nconcat_spikes = [];
        start = 0;
        for i = 1:length(fitblocks)
            block_spikes = spikes(spikes > fitblocks(i) & spikes < fitblocks(i)+fitmovie_seconds_per_block);
            nconcat_spikes = [nconcat_spikes; block_spikes-fitblocks(i)+start];
            start = start + fitmovie_seconds_per_block;
        end
        nspikes{cell} = nconcat_spikes;
        for i = 1:length(testblocks)
            ntspikes{cell}{i} = spikes(spikes > testblocks(i) & spikes < testblocks(i)+testmovie_seconds_per_block) - testblocks(i);
        end
    end

    [STA, center_STA] = STA_Test(concat_spikes, double(fitmovie), true); 
    fittedGLM{i_cell} = glm_fit(concat_spikes, fitmovie, round(center), 'WN_STA', squeeze(STA), 'neighborspikes', 0);
    
    glm_predict(fittedGLM{i_cell}, testmovie, 'testspikes', tspikes, 'neighborspikes', ntspikes)
    
end

% save('/Volumes/Lab/Users/Nora/amacrine_coupled.mat','fittedGLM')
