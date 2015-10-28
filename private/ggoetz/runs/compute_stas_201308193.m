clear;

parpool([1 32])
N_SPIKES_STA = 10000;

% %% data000 - working example
% 
% datarunpath = '/Volumes/Lab/Projects/mining/output/2013-08-19-3/data000-debug2/data000-debug2';
% samples_to_frames = '/Volumes/Lab/Projects/mining/output/2013-08-19-3/data000-debug2/data000-debug2.stf';
% moviechunksfolder = '/Volumes/Lab/Projects/vstim-unpack/unpacked/white-noise/RGB_10_0.48_11111' ;
% interval = 2;
% 
% % STA parameters
% headerCapacity = int32(10000);
% width = int32(32);
% height = int32(64);
% staDepth = int32(76);
% staOffset = int32(0);
% stixelwidth = 10;
% stixelheight = 10;
% refreshtime = 3;
% stafilepath = '/Volumes/Lab/Projects/mining/output/2013-08-19-3/data000-debug2/data000-debug2.sta';
% staFile = edu.ucsc.neurobiology.vision.io.STAFile(stafilepath, headerCapacity, width, height, staDepth, staOffset, stixelwidth, stixelwidth, refreshtime);
% 
% % Load the data run
% datarun = load_data(datarunpath);
% datarun = load_neurons(datarun);
% 
% % RRS works in samples, but it is more natural to work in samples for us.
% % Let's convert the datarun to samples.
% datarun = convert_datarun_times_to_samples(datarun);
% 
% % Get the time of the image refreshes from the ttls
% t_frames = time_imrefresh_from_ttls(datarun.triggers);
% 
% % Progress bar
% fprintf('Calculating STAs.\n');
% fprintf([repmat('.',1,80) '\n']);
% ndots = 0;
% 
% % Get the STAs
% ncells = length(datarun.cell_ids);
% for k = 1:ncells
%     % Update progress bar
%     if (k/ncells)*80 > ndots
%         fprintf('.');
%         ndots = ndots + 1;
%     end
%     
%     % Get spike times and trim
%     st = datarun.spikes{k};
%     if length(st) > N_SPIKES_STA
%         st = st(1:N_SPIKES_STA);
%     end
%     
%     % Get cell ID
%     cellid = datarun.cell_ids(k);
%     
%     % Calculate STA frame indices
%     staind = compute_ta_ind(st, samples_to_frames);
%     
%     % Remap frame indices to movie indices
%     % staindremapped = remap_event_indices(staind, ...
%     %     1:length(t_frames), ...
%     %     reshape(repmat(1:(length(t_frames)/interval), interval, 1), length(t_frames), 1).');
%     % For WN remapping is simply dividing by interval
%     staindremapped = ceil(staind/interval);
%     
%     % Calculate STA from movie indices
%     [sta, e_sta] = compute_ta_from_ind(staindremapped, moviechunksfolder);
% 
%     % Convert the cell array STA to a Vision STA
%     vsta = cell_array_to_vision_sta(sta, e_sta, refreshtime);
%     
%     % Add the STA to the STA file
%     staFile.addSTA(cellid, vsta)
% end
% fprintf('\nSTA calculation done.\n');
% staFile.close()

%% data000 - new test

moviechunksfolder = '/Volumes/Lab/Projects/vstim-unpack/unpacked/white-noise/RGB_10_0.48_11111' ;

% Dataset parameters
datarunpath = '/Volumes/Lab/Projects/mining/output/2013-08-19-3/data000-debug/data000-debug';
interval = 2;

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
% Note reversed width/height compared to what you'd expect.
width = int32(32);
height = int32(64);
staOffset = int32(0);
stixelwidth = 10;
stixelheight = 10;

% Matlab STA time granularity (used by Vision to scale time axis of the STA).
% Corresponds to default STA step size in compute_ta_ind, 120 samples.
% This time granularity should be specified in milliseconds.
refreshtime = 6; 
% Matlab STA depth. Default value is 50 (corresponds to 6000 samples, see
% documentation of compute_ta_ind
staDepth = int32(50);

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
parfor k = 1:ncells
    % Update progress bar - doesn't work with parfor...
    if mod(k, round(ncells/80)) == 0
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
    
    % Convert the cell array STA to a Vision STA
    vsta = cell_array_to_vision_sta(sta, e_sta, refreshtime, stixelwidth);
    
    % Add the STA to the STA file
    staFile.addSTA(cellid, vsta)
end
    
staFile.close()


% %% data000 - working example
% 
% datarunpath = '/Volumes/Lab/Projects/mining/output/2013-08-19-3/data000-debug2/data000-debug1';
% samples_to_frames = '/Volumes/Lab/Projects/mining/output/2013-08-19-3/data000-debug2/data000-debug1.stf';
% moviechunksfolder = '/Volumes/Lab/Projects/vstim-unpack/unpacked/white-noise/RGB_10_0.48_11111' ;
% interval = 2;
% 
% % STA parameters
% headerCapacity = int32(10000);
% width = int32(32);
% height = int32(64);
% staDepth = int32(76);
% staOffset = int32(0);
% stixelwidth = 10;
% stixelheight = 10;
% refreshtime = 3;
% stafilepath = '/Volumes/Lab/Projects/mining/output/2013-08-19-3/data000-debug2/data000-debug2.sta';
% staFile = edu.ucsc.neurobiology.vision.io.STAFile(stafilepath, headerCapacity, width, height, staDepth, staOffset, stixelwidth, stixelwidth, refreshtime);
% 
% % Load the data run
% datarun = load_data(datarunpath);
% datarun = load_neurons(datarun);
% 
% % RRS works in samples, but it is more natural to work in samples for us.
% % Let's convert the datarun to samples.
% datarun = convert_datarun_times_to_samples(datarun);
% 
% % Get the time of the image refreshes from the ttls
% t_frames = time_imrefresh_from_ttls(datarun.triggers);
% 
% % Progress bar
% fprintf('Calculating STAs.\n');
% fprintf([repmat('.',1,80) '\n']);
% ndots = 0;
% 
% % Get the STAs
% ncells = length(datarun.cell_ids);
% for k = 1:ncells
%     % Update progress bar
%     if (k/ncells)*80 > ndots
%         fprintf('.');
%         ndots = ndots + 1;
%     end
%     
%     % Get spike times and trim
%     st = datarun.spikes{k};
%     if length(st) > N_SPIKES_STA
%         st = st(1:N_SPIKES_STA);
%     end
%     
%     % Get cell ID
%     cellid = datarun.cell_ids(k);
%     
%     % Calculate STA frame indices
%     staind = compute_ta_ind(st, samples_to_frames);
%     
%     % Remap frame indices to movie indices
%     % staindremapped = remap_event_indices(staind, ...
%     %     1:length(t_frames), ...
%     %     reshape(repmat(1:(length(t_frames)/interval), interval, 1), length(t_frames), 1).');
%     % For WN remapping is simply dividing by interval
%     staindremapped = ceil(staind/interval);
%     
%     % Calculate STA from movie indices
%     [sta, e_sta] = compute_ta_from_ind(staindremapped, moviechunksfolder);
% 
%     % Convert the cell array STA to a Vision STA
%     vsta = cell_array_to_vision_sta(sta, e_sta, refreshtime);
%     
%     % Add the STA to the STA file
%     staFile.addSTA(cellid, vsta)
% end
% fprintf('\nSTA calculation done.\n');
% staFile.close()

%% data001

%% Calculate the order in which the frames of the natural movie were shown

% repeats = 60;
% frames_a = 3600;
% frames_b = 3600*2;
% fo = frames_order_natural_movie(frames_a, frames_b, repeats);

delete(gcp)