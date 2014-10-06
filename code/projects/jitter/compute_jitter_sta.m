function jitter_stas = compute_jitter_sta(datarun,offset_index,cell_id,bufMov,desired_sta_frames,params)
% COMPUTE_JITTER_STA     compute STA of jittered movie, return a separate STA for each jitter
%
% usage:  jitter_stas = compute_jitter_sta(datarun,offset_index,cell_id,bufMov,desired_sta_frames,params)
%
% arguments:  datarun - datarun struct
%        offset_index - list of offsets
%             cell_id - 
%              bufMov - 
%  desired_sta_frames - 
%              params - struct of optional parameters (see below)
%
% outputs:   jitter_stas - cell array of STA movies
%
%
% optional fields in params, their default values, and what they specify:
%
% spike_fraction    1       fraction of the cell spikes to use
% rand_seed         11111   random seed (used to determine which subset of spikes to use)
%
%


% SET UP OPTIONAL ARGUMENTS

% if not specified, make params empty
if ~exist('params','var');params = [];end

% specify default parameters
defaults.spike_fraction = 1;
defaults.rand_seed = 11111;

% combine user and default parameters
params = default_params( defaults, params);




% BODY OF THE FUNCTION

%% get required info:
% spike times, trigger times, stimulus interval, mdf file, which movie frames have the desired offset


import('edu.salk.snl.cellfinder.cellfinder.*');
import('edu.salk.snl.cellfinder.statistics.*');
import('edu.salk.snl.cellfinder.data.*');

% load neurons file
neuronFile = edu.ucsc.neurobiology.vision.io.NeuronFile(datarun.names.rrs_neurons_path);

% get number of samples
num_samples = neuronFile.getHeader.nSamples;
neuronFile.close;

% get sampling rate
sampling_rate = datarun.sampling_rate;

% get samples per movie frame
samples_per_refresh = NeuronsData.calculateRigFrameTime(datarun.triggers*sampling_rate)/1000*sampling_rate;
samples_per_movie_frame = samples_per_refresh * datarun.stimulus.interval;

% get stimulus parameters
stixel_width = datarun.stimulus.stixel_width;
stixel_height = datarun.stimulus.stixel_height;
stim_interval = datarun.stimulus.interval;
mdf_seed = datarun.stimulus.seed;



%% IDENTIFY THE OFFSET OF EACH FRAME


% compute upper bound on the number of frames
max_num_frames = ceil(length(datarun.triggers)*100/datarun.stimulus.interval);

% initialize variable to store offsets
rand_nums_x = zeros(max_num_frames,1);
rand_nums_y = zeros(max_num_frames,1);

% initialize random number generator
rng = edu.ucsc.neurobiology.vision.math.RandomJavaV2(mdf_seed);

% get seeds
for rr = 1:size(rand_nums_x,1)
    rand_nums_x(rr,1) = rng.nextBits(16);
    rand_nums_y(rr,1) = rng.nextBits(16);
end

% compute offset of each frame
offsets_x = mod(rand_nums_x,stixel_width) - floor(stixel_width/2);
offsets_y = mod(rand_nums_y,stixel_height) - floor(stixel_height/2);




%% LOOP THROUGH OFFSETS AND STA FRAMES TO GET RELEVANT SPIKE TIMES OF EACH CONDITION


% keep track of how many spikes there are in total
max_spike_count = 0;
min_spike_count = Inf;

% seed matlab random number generator if a subset of spikes will be taken
if params.spike_fraction < 1
    rand('twister',params.rand_seed)
end

% loop offsets
for oo = 1:size(offset_index,2)

    % identify frames with the desired offset
    relevant_frames = find((offsets_x==offset_index(1,oo)).*(offsets_y==offset_index(2,oo)));

    % loop STA frames
    for ff = 1:length(desired_sta_frames)

        %% keep only spikes that follow relevant frames by delta t

        % get interval start points (in samples)
        intervals_start = datarun.triggers(1)*sampling_rate + ... time when the first frame started
            (relevant_frames-1)*samples_per_movie_frame + ... plus times when the relevant frames started
            -desired_sta_frames(ff)*samples_per_movie_frame; % plus the delta t offset

        % get interval end times
        intervals_end = intervals_start + samples_per_movie_frame;

        % combine interval start and end times
        interval_edges = sort([intervals_start' intervals_end']);

        % add an edge very early and very late
        interval_edges = [-10*sampling_rate   interval_edges   num_samples+10*sampling_rate];

        % get the spike times (and convert to row vector)
        spikes = datarun.spikes{get_cell_indices(datarun,cell_id)}'*sampling_rate;

        % identify which bin each spike falls in
        [junk,bin] = histc(spikes,interval_edges);

        % get spikes from the even numbered bins
        spikes_to_keep{oo}{ff} = spikes(mod(bin,2)==0);
        
        % take subset of spikes (if desired)
        if params.spike_fraction < 1
            % get the number of spikes in total, and the number to keep
            num_spikes = length(spikes_to_keep{oo}{ff});
            num_to_keep = round(params.spike_fraction*num_spikes);
            % keep only a random subset of them
            rand_ordering = randperm(num_spikes);
            spikes_to_keep{oo}{ff} = sort(spikes_to_keep{oo}{ff}(rand_ordering(1:num_to_keep)));
        end

        % note if this increased the maximum number of spikes
        max_spike_count = max(max_spike_count,length(spikes_to_keep{oo}{ff}));

        % note if this decreased the minimum number of spikes
        min_spike_count = min(min_spike_count,length(spikes_to_keep{oo}{ff}));

    end
end



%% COMPUTE STAS

% put all spike times into a single matrix

% initialize matrix
spikes_matrix = zeros(numel(spikes_to_keep),max_spike_count,'int32');

% counter to note which line is being updated
line_num = 1;

% loop offsets
for oo = 1:size(offset_index,2)
    % loop STA frames
    for ff = 1:length(desired_sta_frames)
        % insert spike times
        spikes_matrix(line_num,1:length(spikes_to_keep{oo}{ff})) = spikes_to_keep{oo}{ff};
        % update line counter
        line_num = line_num + 1;
    end
end

%fprintf('\ncomputing white noise movie...')
% generate white noise movie
%bufMov = ReceptiveField.calculateMovie(mdf_file,int32(datarun.triggers*sampling_rate));


% note what is being computed
fprintf('computing %d STAs for cell id %d (%d to %d spikes per STA)...',...
    size(spikes_matrix,1),cell_id,min_spike_count,max_spike_count)
start_loading = clock; % note when it started

% compute STAs
stas = ReceptiveField.calculateManySTAs(spikes_matrix,int32(datarun.triggers(1)*sampling_rate),bufMov);

% note how long it took
fprintf(' done (%0.1f minutes)\n',etime(clock,start_loading)/60);


%% ASSEMBLE RESULTS

% grab movie parameters: height, width, duration
field_height = double(stas(1).getHeight);
field_width = double(stas(1).getWidth);
num_frames = double(stas(1).getSTADepth);

% counter to note which STA is being loaded
line_num = 1;

% go throuch each offset

% loop offsets
for oo = 1:size(offset_index,2)

    % construct a 4-d matrix of the results
    
    % loop STA frames
    for ff = 1:length(desired_sta_frames)
        
        switch 1 % once method 2 is tested, it should be used instead of method 1
            case 1
                % grab entire buffer of the desired frame
                sta_frame = double(stas(line_num).getSTAFrame(num_frames + desired_sta_frames(ff) - 1).getBuffer);

                % reshape into proper image
                sta_frame = reshape(sta_frame,3,field_width,field_height);
                sta_frame = permute(sta_frame,[3 2 1]);
                
            case 2
                % get desired frame
                sta_frame = sta_from_java_sta(stas(line_num),'frames',num_frames + desired_sta_frames(ff));
                
        end
                
        jitter_stas{oo}(:,:,:,ff) = sta_frame;
        
        % update line counter
        line_num = line_num + 1;
    end
end

