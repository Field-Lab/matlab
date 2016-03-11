clear
Analysis_Path = '/Volumes/Analysis/2015-05-27-3/data001-data005';
datarun_class = load_data([Analysis_Path '/data001/data001'], struct('load_neurons', 0, 'load_params', 1));
model_fit_data = load_data([Analysis_Path '/data002/data002'], struct('load_neurons', 1, 'load_params', 0));

%%
%pd= interleaved_data_prep(model_fit_data, [96000+2300 48000+1700], 5);
block_starts = [1 499 1481 1983 2965 3467 4449 4951 5933 6435];
WN_blocks = block_starts(1:2:end);
NSEM_block_starts = block_starts(2:2:end);

%%
% block frames is a vector of the different block lengths,
% ex = [1200*3 1200] aka even fitting blocks first and odd testing blocks
% second. Weird order is so that frames(block i) = block_frames(mod(i,2)+1)
%
% Output:
% preppeddata structure, including start_time which has the starts of the
% blocks and p which are the optional parameters you included
% It can also include: fitspikes, testspikes, fitmovie, testmovie, etc
% which depend on which options you choose. They will be in cells and each
% cell corresponds to a block of data.
%
% cell_spec
% Setting cell_specification will return the organized spikes for those
% cells. Details for what cell_specification should be are detailed in the
% get_cell_indices function

%% Find Block Times

cell_spec = [202];
visual_check = 1;

for i_trigger = WN_blocks(3)
    % load up the triggers and find triggers per block
    triggers = model_fit_data.triggers(i_trigger+(1:1000));
    block_frames = [3600 1200];
    n_repeats = 10;
    
    % calculate some things
    repeat_only = false;
    n_blocks = n_repeats*2;
    testmovie_only = 0;
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
    if visual_check
        figure;
        plot(diff(triggers))
        hold on
        plot(start_trigger, 1/1.2*ones(length(start_trigger)), '*')
        title('Identifying the Block Starts')
    end
    
    % find the start times
    start_time = triggers(start_trigger);
    prepped_data.start_time = start_time;
    
    % Find and organize the spikes
    disp(['Trigger finding took ' num2str(toc) ' seconds.'])
    
    % Organize the cell spikes
    if cell_spec
        disp('Loading spikes...')
        tic
        if isstruct(datarun_class)
            cids = get_cell_indices(datarun_class, cell_spec);
        else
            cids = get_cell_indices(model_fit_data, cell_spec);
        end
        spikes = cell(n_blocks,length(cids)); % initialize spike cell array
        for i_cell = 1:length(cids)
            spikes_temp = model_fit_data.spikes{cids(i_cell)};
            for i_trigger = 1:n_blocks
                block_length = block_frames(mod(i_trigger,2)+1)/120;
                spikes{i_trigger, i_cell} = spikes_temp( (spikes_temp > start_time(i_trigger)) & (spikes_temp < (start_time(i_trigger) + block_length)) ) - start_time(i_trigger);
            end
        end
        if repeat_only
            prepped_data.testspikes = spikes;
        else
            prepped_data.fitspikes = spikes(2:2:end, :);
            prepped_data.testspikes = spikes(1:2:end, :);
        end
        clear spikes
        
        % Visual check
        if visual_check
            figure;
            hold on
            for i = 1:size(prepped_data.testspikes, 1)
                plot(prepped_data.testspikes{i, 1}, i*ones(length(prepped_data.testspikes{i, 1})), 'k.')
            end
            title('Checking Cell Repeats')
        end
        disp(['Spike organization took ' num2str(toc) ' seconds.'])
    end
    
end

%% Load up and organize the stimulus

%
% stimulus_name
% can be an XML type specification (ie BW-8-1-0.48-11111) currently only
% works for BW, NOT RGB
% Or the filename of a rawmovie
%
% seed
% the fixed seed for white noise repeats
%
% variable seed
% the original variable seed for changing white noise segments. Usually the
% same as fixed seed.
%
% STA_check
% use the loaded stimulus and the spikes to calculate the STA of the first
% cell
%
% RGN_params_A_M_C
% Shouldn't need to change these, but they are the parameters of the random
% nubmer generator for changing the variable seed
%
% activity_movie
% create a movie with the stimulus overlaid with the cell activity
%
% testmovie only
% only load up the test movie (much faster)


% loading up the appropriate stimulus
if p.Results.stimulus_name
    disp('Loading stimulus...')
    tic
    
    % If it's an XML movie
    if strcmp(p.Results.stimulus_name(1:2), 'BW')
        seed_fixed = p.Results.seed;
        
        % Testmovie is easy
        prepped_data.testmovie = get_WN_movie(['/Volumes/Analysis/stimuli/white-noise-xml/' p.Results.stimulus_name '-0.48-' num2str(seed_fixed) '.xml'], block_frames(2));
        
        if ~testmovie_only
            % RNG parameters for updating seed for novel blocks. Shouldn't
            % really have to change these from defaults very often if at
            % all
            a = p.Results.RGN_params_A_M_C(1);
            m = p.Results.RGN_params_A_M_C(2);
            c = p.Results.RGN_params_A_M_C(3);
            
            % Load up each novel block
            seed_variable = p.Results.variable_seed_start;
            prepped_data.fitmovie = cell(floor(n_blocks/2),1);
            for i = 1:floor(n_blocks/2)
                seed_variable = mod( (a*seed_variable + c), m);
                prepped_data.fitmovie{i} = get_WN_movie(['/Volumes/Lab/Users/akheitman/NSEM_Home/Stimuli/BW-8-1-0.48-11111_RNG_16807/xml_files/novel/' p.Results.stimulus_name '-0.48-' num2str(seed_variable) '.xml'], block_frames(1));
            end
        end
        
        % Need to set up RGB
    elseif strcmp(p.Results.stimulus_name(1:2), 'RG')
        error('this code is not set up for RGB movies yet')
        
        % NSbrownian and LPF movies
    elseif strcmp(p.Results.stimulus_name(1:2), 'NS') || strcmp(p.Results.stimulus_name(1:10),'lpf')
        if strcmp(p.Results.stimulus_name(1:10),'NSbrownian')
            movie_path = '/Volumes/Data/Stimuli/movies/eye-movement/current_movies/NSbrownian/matfiles/';
            if strcmp(p.Results.stimulus_name(12), '3')
                switch p.Results.stimulus_name(17);
                    case 'A'
                        offset = 0;
                    case 'B'
                        offset = 3000;
                    case 'C'
                        offset = 6000;
                    case 'D'
                        offset = 9000;
                    case 'E'
                        offset = 12000;
                    case 'F'
                        offset = 15000;
                end
            elseif strcmp(A(12), '6')
                offset = 0;
            else
                error('Not set up for this movie')
            end
        elseif strcmp(p.Results.stimulus_name(1:10),'lpf')
            movie_path = '/Volumes/Data/Stimuli/movies/lpf/version-2/njb_code_and_movies/matfiles/';
            offset = 0;
        else
            error('Not set up for this movie')
        end
        % testmovie
        block_seconds = block_frames/120;
        prepped_data.testmovie = zeros(160, 320, block_frames(2), 'uint8');
        idx = 1:120;
        for i_chunk = (1:block_seconds(2))+offset
            chunk = load([movie_path 'movie_chunk_' num2str(i_chunk) '.mat']);
            prepped_data.testmovie(:,:,idx) = permute(chunk.movie, [2, 1, 3]);
            idx = idx+120;
        end
        if ~testmovie_only
            % fitmovie
            offset = block_seconds(2);
            for i_trigger = 1:floor(n_blocks/2)
                prepped_data.fitmovie{i_trigger} = zeros(160, 320, block_frames(1), 'uint8');
                idx = 1:120;
                for i_chunk = (1:block_seconds(2))+offset
                    chunk = load([movie_path 'movie_chunk_' num2str(i_chunk) '.mat']);
                    prepped_data.fitmovie{i_trigger}(:,:,idx) = permute(chunk.movie, [2, 1, 3]);
                    idx = idx+120;
                end
                offset = offset + block_seconds(1);
            end
            clear chunk
        end
        
        % Arbitrary other rawMovie
    elseif strcmp(p.Results.stimulus_name(end-8:end), '.rawMovie')
        % get beginning of movie for the raster
        prepped_data.testmovie = get_rawmovie(p.Results.stimulus_name, block_frames(2), 0);
        offset = block_frames(2);
        % get the rest of the movie in blocks
        if ~testmovie_only
            for i_trigger = 1:1:floor(n_blocks/2)
                prepped_data.fitmovie{i_trigger} = get_rawmovie(p.Results.stimulus_name(end-8:end), block_frames(1), offset);
                offset = offset + block_frames(1);
            end
        end
    else
        error('Not set up for this movie');
    end
    disp(['Loading the stimulus took ' num2str(toc) ' seconds.'])
end

% Do an STA check using some spikes and whatnot
if p.Results.STA_check && ~testmovie_only
    disp('Doing an STA check...')
    tic
    if ~iscell(prepped_data.fitspikes) || ~iscell(prepped_data.fitmovie)
        warn('Must load stimulus and spikes to do an STA check')
    else
        prepped_data.STA = STA_from_blocks(prepped_data.fitspikes,prepped_data.fitmovie);
    end
    disp(['STA took '  num2str(toc) ' seconds.'])
end

if ischar(p.Results.activity_movie)
    disp('Making the activity movie...')
    tic
    if ~ischar(p.Results.cell_spec)
        warning('You did not specify a cell type. The movie will just have the cell ids you listed')
    end
    
    for i_cell = 1:length(cids)
        spikes_frame = floor(cell2mat(prepped_data.testspikes(:,i_cell)) * 120);
        for i_frame = 1:block_frames(2)
            res.spikes(i_cell, i_frame) = sum(spikes_frame == i_frame);
        end
        res.centers(i_cell,:) = p.Results.datarun_class.vision.sta_fits{cids(i_cell)}.mean;
    end
    figure;
    res_spikes_plot(prepped_data.testmovie, res, p.Results.activity_movie, 'scaling', 4)
    disp(['Making the movie took ' num2str(toc) '  seconds.'])
end
