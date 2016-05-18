%% NSEM fitting and testing only for RPE piece 2015-05-27-11

%% Load data
datapath = '2015-05-27-11/data001-data005-norefit/data002-from-data001_data002_data003_data004_data005/data002-from-data001_data002_data003_data004_data005';
datarun= load_data(datapath);
datarun = load_neurons(datarun);
monitor_refresh = 100/median(diff(datarun.triggers));

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

clear WN4 NSEM repeats_within_block repeat_starts block_starts

% On par
cells = get_cell_indices(datarun_class, 'On Parasol');

%% Load up movies

disp('Loading stimulus...')
tic

% If it's an XML movie
seed_fixed = 11111;
block_frames = [1200 3600];

% Testmovie is easy
prepped_data.testmovie = get_WN_movie(['/Volumes/Analysis/stimuli/white-noise-xml/BW-8-1-0.48-' num2str(seed_fixed) '.xml'], block_frames(2));

% RNG parameters for updating seed for novel blocks. Shouldn't
% really have to change these from defaults very often if at
% all
a = 16807;
m =  2^31 - 1;
c = 0;

% Load up each novel block
seed_variable = seed_fixed;
to_use = zeros(10,10);
to_use(:,2:2:end) = 1;
to_use = to_use(:);
count = 0;
for i = 1:length(to_use)
    seed_variable = mod( (a*seed_variable + c), m);
    if to_use(i)
        count = count+1;
        try
            prepped_data.fitmovie{count} = get_WN_movie(['/Volumes/Lab/Users/Nora/Stimulus/BW_XML/BW-8-1-0.48-' num2str(seed_variable) '.xml'], block_frames(1));
        catch
            eval(['! /Volumes/Lab/Users/Nora/Stimulus/movie-xml-maker BW-8-1-0.48-' num2str(seed_variable)]) % OS system command to make the XML file
            prepped_data.fitmovie{count} = get_WN_movie(['/Volumes/Lab/Users/Nora/Stimulus/BW_XML/BW-8-1-0.48-' num2str(seed_variable) '.xml'], block_frames(1));
        end
    end
end


%% Load up cell info

for i_cell = 1%:length(cells)
    
    disp(i_cell)
    cell = cells(i_cell);
    spikes = datarun.spikes{cell};
    % Concatenate fit spikes
    fitblocks = WN8(2:2:end);
    fitmovie_frames_per_block = 3600;
    fitmovie_seconds_per_block = fitmovie_frames_per_block/monitor_refresh;
    concat_spikes = [];
    start = 0;
    for i = 1:length(fitblocks)
        block_spikes = spikes(spikes > fitblocks(i) & spikes < fitblocks(i)+fitmovie_seconds_per_block);
        concat_spikes = [concat_spikes; block_spikes-fitblocks(i)+start];
        start = start + fitmovie_seconds_per_block;
    end
    % Organize test spikes
    testblocks = WN8(1:2:end);
    testmovie_frames_per_block = 1200;
    testmovie_seconds_per_block = testmovie_frames_per_block/monitor_refresh;
    for i = 1:length(testblocks)
        tspikes{i} = spikes(spikes > testblocks(i) & spikes < testblocks(i)+testmovie_seconds_per_block) - testblocks(i);
    end
    
    center = datarun_class.vision.sta_fits{cell}.mean;
    center(2) = 40 - center(2);
    
    %{
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
    %}

    %[STA, center_STA] = STA_Test(concat_spikes, double(fitmovie), true); 
    fittedGLM{i_cell} = glm_fit(concat_spikes, fitmovie, round(center), 'neighborspikes', 0);
    
    glm_predict(fittedGLM{i_cell}, testmovie, 'testspikes', tspikes, 'neighborspikes', ntspikes)
    
end