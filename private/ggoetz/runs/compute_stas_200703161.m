clear;

parpool([1 32])
addpath(genpath('/home/ggoetz/Research/code/common-chichilnisky-lab/matlab/private/ggoetz'));
addpath(genpath('/home/ggoetz/Research/code/common-chichilnisky-lab/matlab/utilities'));
N_SPIKES_STA = 5000;
zeroval = 128; % For WN movies, this should be 0.5 and for raw movies, 128

%% data003

moviechunksfolder = '/Volumes/Lab/Projects/vstim-unpack/unpacked/np/npg-128-64-64-16-[-0_5]-[-1_0]';
% % Optional: do it once to convert a raw movie to mat chunks.
% % Once you've done for a movie, no need to convert to chunks ever again,
% % so uncomment the following three lines.
% moviepath = '/Volumes/Data/Stimuli/movies/np/npg-128-64-64-16-[-0_5]-[-1_0]';
% greyscale = true;
% unpack_rawmovie(moviepath, moviechunksfolder, greyscale);

% Dataset parameters
datarunpath = '2007-03-16-1/data003';
interval = 4;

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

% % Use that to cache a map of samples to frame indices. 
% % This only needs to be calculated once per dataset and it's slow, so 
% % comment out the following two lines if you need to run the sta 
% % calculation more than once.
% map_samples_to_frames(1:length(t_frames), t_frames, datarun.duration, samples_to_frames);

% Vision STA parameters
headerCapacity = int32(10000);
width = int32(128);
height = int32(64);
staOffset = int32(0);
stixelwidth = 5;
stixelheight = 5;

% Matlab STA time granularity (used by Vision to scale time axis of the STA).
% Corresponds to default STA step size in compute_ta_ind, 120 samples.
% This time granularity should be specified in milliseconds.
refreshtime = 6; 
% Matlab STA depth. Default value is 50 (corresponds to 6000 samples, see
% documentation of compute_ta_ind
staDepth = int32(50);

% Instantiate Vision STA file
staFile = edu.ucsc.neurobiology.vision.io.STAFile(stafilepath, headerCapacity, width, height, staDepth, staOffset, stixelwidth, stixelwidth, refreshtime);

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
parfor k = 1:ncells
    % Update progress bar - doesn't work with parfor...
    if mod(length(dir(stastempfolder)) - 2, round(ncells/80)) == 0
        fprintf('.');
    end
    
    % Get spike times and trim
    st = datarun.spikes{k};
    if length(st) > N_SPIKES_STA
        st = st(1:N_SPIKES_STA);
    end
    
    % Get cell ID
    cellid = datarun.cell_ids(k);
    
    % Calculate STA frame indices
    staind = compute_ta_ind(st, samples_to_frames);
    
    % Remap frame indices to movie indices
    % This is where the interval matters.
    % You can also remap more complex visual stimuli like natural movies
    % here.
    staindremapped = ceil(staind/interval);
    
    % Calculate STA from movie indices
    [sta, e_sta] = compute_ta_from_ind(staindremapped, moviechunksfolder, zeroval);
    
    % Save the result
    save_parfor_stas(fullfile(stastempfolder, sprintf('sta_%s.mat', num2str(cellid))), ...
        sta, e_sta);
end
fprintf('\nSTA calculation done. Saving...\n');

for k = 1:20
    cellid = datarun.cell_ids(k);
    load(fullfile(stastempfolder, sprintf('sta_%s.mat', num2str(cellid))))
    
    % Convert the cell array STA to a Vision STA
    vsta = cell_array_to_vision_sta(sta, e_sta, refreshtime, stixelwidth);
    
    % Add the STA to the STA file
    staFile.addSTA(cellid, vsta)
end
    
staFile.close()
delete(gcp);
