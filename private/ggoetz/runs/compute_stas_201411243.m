clear;
addpath(genpath('/home/ggoetz/Research/code/common-chichilnisky-lab/matlab/private/ggoetz'));
addpath(genpath('/home/ggoetz/Research/code/common-chichilnisky-lab/matlab/utilities'));

%% Unpack raw movies if necessary
% Optional: do it once to convert a raw movie to mat chunks.
% Once you've done for a movie, no need to convert to chunks ever again,
% so uncomment the following three lines.

% moviechunksfolder = '/Volumes/Lab/Projects/vstim-unpack/unpacked/eye-movement/NSbrownian_3000_A_025';
% moviepath = '/Volumes/Data/Stimuli/movies/eye-movement/current_movies/NSbrownian/NSbrownian_3000_movies/NSbrownian_3000_A_025.rawMovie';
% greyscale = true;

% unpack_rawmovie(moviepath, moviechunksfolder, greyscale);


%% STA calculations setup

% parpool([1 32])
N_SPIKES_STA = 15000;

%% data009

moviechunksfolder = '/Volumes/Lab/Projects/vstim-unpack/unpacked/eye-movement/NSbrownian_3000_A_025' ;

% Dataset parameters
datarunpath = '/Volumes/Analysis/2014-11-24-3/mv/data009/data009';
interval = 1;
repeats = 50;
frames_a = 3600;
frames_b = 3600*2;
wrap = 3000*120;

% Load datarun
datarun = load_data(datarunpath);
datarun = load_neurons(datarun);

% Set paths
samples_to_frames = [datarun.names.rrs_prefix '.stf'];
stafilepath = [datarun.names.rrs_prefix '.sta'];

% RRS works in samples, but it is more natural to work in samples for us.
% Let's convert the datarun to samples.
datarun = convert_datarun_times_to_samples(datarun);

% Get the time of the image refreshes from the ttls
t_frames = time_imrefresh_from_ttls(datarun.triggers);

% Calculate the order in which the natural movie frames were shown
fo = frames_order_natural_movie(frames_a, frames_b, repeats, wrap);
% Stimulus appears to have been cut short, so we're going to trim the frame
% order so that things match
fo = fo(1:length(t_frames));

% Use that to cache a map of samples to frame indices. 
% This only needs to be calculated once per dataset and it's slow, so 
% comment out the following two lines if you need to run the sta 
% calculation more than once.
% map_samples_to_frames(fo, t_frames, datarun.duration, samples_to_frames);

% Calculate when the repeats occured
tr = get_start_end_time_repeats(fo, t_frames);

% Vision STA parameters
headerCapacity = int32(10000);
% Note reversed width/height compared to what you'd expect.
width = int32(160);
height = int32(320);
staOffset = int32(0);
stixelwidth = 2;
stixelheight = 2;

% Matlab STA time granularity 
matlab_sta_times = [6000, 1200, 240];  % Specify here if you want some other value that the default (small) time step
% Granularity in milliseconds for Vision to scale the time axis correctly.
refreshtime = 12; 
% Matlab STA depth. Default value is 50 (corresponds to 6000 samples, see
% documentation of compute_ta_ind
staDepth = int32(30);

% Instantiate Vision STA file
staFile = edu.ucsc.neurobiology.vision.io.STAFile(stafilepath, headerCapacity, width, height, staDepth, staOffset, stixelwidth, stixelheight, refreshtime);

% STA temp folder - needed to work around clunkiness of Matlab parallel
% computations.
stastempfolder = split(datarun.names.rrs_prefix, filesep);
stastempfolder = fullfile(join(stastempfolder(1:(end-1)), filesep), 'stastemp');
if exist(stastempfolder, 'dir') == 0
    mkdir(stastempfolder);
end

% Progress bar
fprintf('Calculating STAs.\n');
fprintf([repmat('.',1,80) '\n']);
ndots = 0;

% Get the STAs
ncells = length(datarun.cell_ids);
cells_to_delete = zeros(ncells, 1);
for k = 1:ncells
% parfor k = 1:ncells
    % Update progress bar - doesn't work with parfor...
    if mod(k, round(ncells/80)) == 0
        fprintf('.');
    end
    
    % Get spike times
    st = datarun.spikes{k};
    
    % For natural scenes, remove raster spikes from the train
    st = remove_repeats_from_st(tr, st);
    
    % Trim train to the max allowed size
    if length(st) > N_SPIKES_STA
        st = st(1:N_SPIKES_STA);
    end
    
    % Get cell ID
    cellid = datarun.cell_ids(k);
    
    % Calculate STA frame indices
    staind = compute_ta_ind(st, samples_to_frames, matlab_sta_times);
    if isempty(staind)
        cells_to_delete(k) = 1;
        continue;
    end
    
    % Remap frame indices to movie indices
    % This is where the interval matters.
    % You can also remap more complex visual stimuli like natural movies
    % here, using for example remap_event_indices(staing, i_pre, i_post).
    % To get null of your STA distribution: randomly permutate indices
    % here.
    staindremapped = ceil(staind/interval);
    
    % Calculate STA from movie indices
    [sta, e_sta] = compute_ta_from_ind(staindremapped, moviechunksfolder);
    
    % Save the result
    save_parfor_stas(fullfile(stastempfolder, sprintf('sta_%s.mat', num2str(cellid))), ...
        sta, e_sta);
end
fprintf('\nSTA calculation done. Saving...\n');

datarun.cell_ids(logical(cells_to_delete)) = [];
ncells = length(datarun.cell_ids);
for k = 1:ncells
    cellid = datarun.cell_ids(k);
    load(fullfile(stastempfolder, sprintf('sta_%s.mat', num2str(cellid))))
    
    % Possibly try to take out some offset in each sta
    sta = remove_offset_from_sta(sta, mean(sta{1}(:)));
    
    % Convert the cell array STA to a Vision STA
    vsta = cell_array_to_vision_sta(sta, e_sta, refreshtime, stixelwidth);
    
    % Add the STA to the STA file
    staFile.addSTA(cellid, vsta)
end
    
staFile.close()


%% data001



delete(gcp)