clear;

addpath(genpath('/home/ggoetz/Research/code/common-chichilnisky-lab/matlab/private/ggoetz'));
N_SPIKES_STA = 10000;

%% data000 

datarunpath = '/Volumes/Lab/Projects/mining/output/2013-08-19-3/data000-debug2/data000-debug2';
samples_to_frames = '/Volumes/Lab/Projects/mining/output/2013-08-19-3/data000-debug2/data000-debug2.stf';
moviechunksfolder = '/Volumes/Lab/Projects/vstim-unpack/unpacked/white-noise/RGB_10_0.48_11111' ;
interval = 2;

% STA parameters
headerCapacity = int32(10000);
width = int32(32);
height = int32(64);
staDepth = int32(76);
staOffset = int32(0);
stixelwidth = 10;
stixelheight = 10;
refreshtime = 2;
stafilepath = '/Volumes/Lab/Projects/mining/output/2013-08-19-3/data000-debug2/data000-debug2.sta';
staFile = edu.ucsc.neurobiology.vision.io.STAFile(stafilepath, headerCapacity, width, height, staDepth, staOffset, stixelwidth, stixelwidth, refreshtime);

% Load the data run
datarun = load_data(datarunpath);
datarun = load_neurons(datarun);

% RRS works in samples, but it is more natural to work in samples for us.
% Let's convert the datarun to samples.
datarun = convert_datarun_times_to_samples(datarun);

% Get the time of the image refreshes from the ttls
t_frames = time_imrefresh_from_ttls(datarun.triggers);

% Progress bar
fprintf('Calculating STAs.\n');
fprintf([repmat('.',1,80) '\n']);
ndots = 0;

% Get the STAs
ncells = length(datarun.cell_ids);
ncells = 1;
for k = 1:ncells
    % Update progress bar
    if (k/ncells)*80 > ndots
        fprintf('.');
        ndots = ndots + 1;
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
    % staindremapped = remap_event_indices(staind, ...
    %     1:length(t_frames), ...
    %     reshape(repmat(1:(length(t_frames)/interval), interval, 1), length(t_frames), 1).');
    % For WN remapping is simply dividing by interval
    staindremapped = ceil(staind/interval);
    
    % Calculate STA from movie indices
    [sta, e_sta] = compute_ta_from_ind(staindremapped, moviechunksfolder);

    % Convert the cell array STA to a Vision STA
    vsta = cell_array_to_vision_sta(sta, e_sta, refreshtime);
    
    % Add the STA to the STA file
    staFile.addSTA(cellid, vsta)
end
fprintf('\nSTA calculation done.\n');
staFile.close()

%% data001

%% Calculate the order in which the frames of the natural movie were shown

% repeats = 60;
% frames_a = 3600;
% frames_b = 3600*2;
% fo = frames_order_natural_movie(frames_a, frames_b, repeats);