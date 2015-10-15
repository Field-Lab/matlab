function  [start_time, testmovie, fitmovie, testspikes, fitspikes] = interleaved_data_prep(datarun, block_frames, n_repeats, varargin)
% This assumes that the odd blocks are repeated rasters and the even blocks
% are for fitting
%
% block frames is a vector of the different block lengths, in order
% = [1200*3 1200] aka even, fitting blocks first and odd, testing blocks
% second. Weird order is so that frames(block i) = block_frames(mod(i,2)+1)
%
% Optional:
%
% Setting cell_specification will return the organized spikes for those
% cells. Details for what cell_specification should be are detailed in the
% get_cell_indices function
%
% Stimulus name:
% can be the beginning of an XML type specification (ie BW-8-1)
% Or the filename of a rawmovie 
% Or the path to a folder of mat files

p = inputParser;
p.addParameter('cell_spec', 0)
p.addParameter('stimulus_name', 0)
p.addParameter('seed', 11111)
p.addParameter('variable_seed_start', 11111)
p.addParameter('visual_check', 1)
p.addParameter('STA_check', 0)
p.addParameter('RGN_params_A_M_C', [16807, 2^31 - 1, 0]); % Shouldn't really ever have to change this
p.parse(varargin{:});

% Just set to 0 if we aren't using
if ~p.Results.stimulus_name; testmovie = 0; fitmovie = 0; end
if ~p.Results.cell_spec; testspikes = 0; fitspikes = 0; end


%% Find Block Times

% load up the triggers and find triggers per block
triggers = datarun.triggers;

%
if length(block_frames) == 1
    repeat_only = true;
    n_blocks = n_repeats;
    block_frames = [block_frames, block_frames]; % not pretty, but easy
elseif length(block_frames) == 2
    repeat_only = false;
    n_blocks = n_repeats*2;
else
    error('Block_frames input is the wrong length')
end
triggers_per_block = block_frames/100;

% Initialize and define things
start_trigger = zeros(n_blocks,1);
start_trigger(1) = 1;

% Find the trigger numbers of the block starts
for i = 1:n_blocks
    end_trigger = start_trigger(i) +  triggers_per_block(mod(i,2)+1);
    start_trigger(i+1) = end_trigger + 1;
end
start_trigger = start_trigger(start_trigger <= length(triggers));
n_blocks = length(start_trigger);

% visual check
if p.Results.visual_check
    figure;
    plot(diff(triggers))
    hold on
    plot(start_trigger, 1/1.2*ones(length(start_trigger)), '*')
    title('Identifying the Block Starts')
end

% find the start times
start_time = triggers(start_trigger);


%% Find and organize the spikes

% Organize the cell spikes
if p.Results.cell_spec
    cids = get_cell_indices(datarun, p.Results.cell_spec);
    spikes = cell(n_blocks,length(cids)); % initialize spike cell array
    for i_cell = 1:length(cids)
        spikes_temp = datarun.spikes{cids(i_cell)};
        for i_block = 1:n_blocks
            block_length = block_frames(mod(i_block,2)+1)/120;
            spikes{i_block, i_cell} = spikes_temp( (spikes_temp > start_time(i_block)) & (spikes_temp < (start_time(i_block) + block_length)) ) - start_time(i_block);
        end
    end
    if repeat_only
        testspikes = spikes;
        fitspikes = 0;
    else
        fitspikes = spikes(2:2:end, :);
        testspikes = spikes(1:2:end, :);
    end
    clear spikes
    
    % Visual check
    if p.Results.visual_check
        figure;
        hold on
        for i = 1:length(testspikes)
            plot(testspikes{i, 1}, i*ones(length(testspikes{i, 1})), 'k.')
        end
        title('Checking Cell Repeats')
    end
end

%% Load up and organize the stimulus

% loading up the appropriate stimulus
if p.Results.stimulus_name
    
    % If it's an XML movie
    if strcmp(p.Results.stimulus_name(1:2), 'BW')
        seed_fixed = p.Results.seed;
        
        % Testmovie is easy
        testmovie = subR_get_movie(['/Volumes/Analysis/stimuli/white-noise-xml/' p.Results.stimulus_name '-0.48-' num2str(seed_fixed) '.xml'], block_frames(2));
        
        if ~repeat_only
            % RNG parameters for updating seed for novel blocks. Shouldn't
            % really have to change these from defaults very often if at
            % all
            a = p.Results.RGN_params_A_M_C(1);
            m = p.Results.RGN_params_A_M_C(2);
            c = p.Results.RGN_params_A_M_C(3);
            
            % Load up each novel block
            seed_variable = p.Results.variable_seed_start;
            fitmovie = cell(floor(n_blocks/2),1);
            for i = 1:floor(n_blocks/2)
                seed_variable = mod( (a*seed_variable + c), m);
                fitmovie{i} = subR_get_movie(['/Volumes/Lab/Users/akheitman/NSEM_Home/Stimuli/BW-8-1-0.48-11111_RNG_16807/xml_files/novel/' p.Results.stimulus_name '-0.48-' num2str(seed_variable) '.xml'], block_frames(1));
            end
        end
        
        % Need to set up NSEM and other movies
    elseif strcmp(p.Results.stimulus_name(1:2), 'RG')
        error('this code is not set up for RGB movies yet')
    elseif strcmp(p.Results.stimulus_name(1:2), 'NS')
        error('this code is not set up for NSEM movies yet')
    end
end

% Do an STA check using some spikes and whatnot
if p.Results.STA_check
    if ~iscell(spikes) || ~iscell(fitmovie)
        warn('Must load stimulus and spikes to do an STA check')
    else
        stim_size = size(fitmovie{1});
        STA = zeros(stim_size(1), stim_size(2), 30); % 30 frame STA
        for i_block = 1:length(fitmovie)
            block_spikes = spikes{2*i_block, 1};
            for i_spike = 1:length(block_spikes)
                spike_frame = floor(120*block_spikes(i_spike));
                if spike_frame > 29
                   STA = STA + fitmovie{i_block}(:,:,(spike_frame-29):spike_frame);
                end
            end
        end
        figure;
        for i = 1:30
           imagesc(STA(:,:,i))
           axis image
           colormap gray
           pause(0.1)
        end
        imagesc(sum(STA(:,:,24:27), 3))
        axis image
        colormap gray
    end
end



end

function movie = subR_get_movie(mdf_file, frames)
% GET_MOVIE      Load a white noise movie

% Make fake triggers for loading
n_triggers = frames/100 + 100; % one extra for good luck
triggers = 0:(1/1.2):((1/1.2)*n_triggers);

% load movie
[mvi] = load_movie(mdf_file, triggers);

% grab movie parameters: height, width, duration
height = double(mvi.getHeight);
width = double(mvi.getWidth);
duration = double(mvi.size);
% refresh = double(mvi.getRefreshTime);

if exist('frames','var')
    if frames<=duration
        tduration=frames;
    else
        error('movie too short');
    end
end
  
movie=zeros(height,width,tduration);
for i=1:tduration
    % grab movie frame
    F = mvi.getFrame(i-1).getBuffer;

    % reshape into proper image
    F = reshape(F,3,width,height);
    F = permute(F,[3 2 1]);
    movie(:,:,i) = F(:,:,1); % just take one channel since BW
end
end


% Extra code
% define the number of blocks and their order. 
% Usually n_block_types will be 2: fitting and testing
% n_block_types = length(block_frames);
% n_blocks = 100; % how should this be figured out?
% if p.Results.block_order
%    block_order = p.Results.block_order;
% else
%   block_order =  repmat(1:2,1, ceil(n_blocks/n_block_types));
% end
