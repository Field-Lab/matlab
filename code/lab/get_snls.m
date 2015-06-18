function datarun = get_snls(datarun, cell_spec, varargin)
% GET_SNLS       compute SNLs and store in datarun, optimized version
% usage:  datarun = get_snls(datarun, cell_spec, <params>)
%
% arguments:  datarun - datarun struct
%           cell_spec - which cells (see get_cell_indices for options)
%            <params> - struct or list of optional parameters (see below)
%
% outputs:    datarun - datarun struct with results stored in datarun.stas.snls{}
%
%
% optional parameters, their default values, and what they specify:
%
%
% verbose           true            show output
% movie             []              java movie object specifying stimulus
%                                       if empty, will look in datarun.stimulus.java_movie
%                                       note: if the movie is for a very high dimensional stimulus,
%                                       it is more efficient to compute the
%                                       movie in advance with
%                                       LOAD_JAVA_MOVIE
%                                       
% frames            ':'             which frames to use
%                                       see parse_frame_spec for options
% new               false           if SNL already exists, compute a new one?
% stimuli           []              how many stimulus frames to use
%                                       if empty, go until end
% 
% start_time        0               time of the first stimulus frame to use (seconds)
% end_time          datarun.duration  time when the stimulus ends
%
%
% 
%
% See also: LOAD_JAVA_MOVIE
%
% 2009-09  gauthier
% 2013-06  phli, making more modular, optimized by changing STRF format to
%                use sparse matrix operations efficiently and to match the
%                vision movie frame buffer format so that the movie frames
%                no longer need to be permuted.  Old code currently
%                retained as GET_SNLS_OLD
%
% 2014-01 gdf, the start_stim and end_stim functionality has been removed
% 
%
% See also: GET_SNLS_OLD, COMPUTE_GEN_SIGNALS
%


% SET UP OPTIONAL ARGUMENTS

p = inputParser;

% specify list of optional parameters
p.addParamValue('verbose', true);
p.addParamValue('movie', []);
p.addParamValue('frames', ':');
p.addParamValue('new', false);
p.addParamValue('stimuli', []);
p.addParamValue('start_time', 0, @(x)numel(x)==1 && x>=0 && x <= datarun.duration);
p.addParamValue('end_time', datarun.duration, @(x)numel(x) == 1);
p.addParamValue('marks', 'marks');
p.addParamValue('simplemarkthresh', 4);

% resolve user input and default values
p.parse(varargin{:});

% get params struct
params = p.Results;




% identify how many cells need SNLs computed

% get cell indices
cell_indices_ = get_cell_indices(datarun,cell_spec);

% ensure proper fields exist in datarun
if ~isfield(datarun.stas,'snls')
    datarun.stas.snls = cell(length(datarun.cell_ids),1);
else
    if length(datarun.stas.snls) < length(datarun.cell_ids)
        datarun.stas.snls{length(datarun.cell_ids)} =[];
    end
end

% if new computation
if params.new
    % compute for all cells
    cell_indices = cell_indices_;
else
    % otherwise, identify only the cells that need it
    cell_indices = [];
    
    % check each cell
    for cc = 1:length(cell_indices_)
        if isempty(datarun.stas.snls{cell_indices_(cc)})
            cell_indices = [cell_indices cell_indices_(cc)]; %#ok<AGROW>
        end
    end
end



% if none, quit
if isempty(cell_indices)
    return
end


% prepare movie

if isempty(params.movie)
    % if none was specified, get movie from datarun
    if isfield(datarun.stimulus,'java_movie') && ~isempty(datarun.stimulus.java_movie)
        movie = datarun.stimulus.java_movie;
    else
        % if datarun doesn't have a movie, give an error
        error('Stimulus movie not specified.  Use datarun = load_java_movie(datarun,<movie_xml_path>);')
    end
else
    % otherwise, use what was provided
    movie = params.movie;
end



 

% get STRFs and spike times

% identify refresh time
refresh_time = movie.getRefreshTime/1000;

% note stimulus size
field_width = movie.getWidth;
field_height = movie.getHeight;

% check for consistency
if ~all([field_width field_height] == [datarun.stimulus.field_width datarun.stimulus.field_height])
    error('')
end

% initialize storage variables
spikes = sparse(zeros(movie.size,length(cell_indices)));
rois = sparse(false(field_width*field_height,length(cell_indices)));
strfs = cell(length(cell_indices),1);

% show output
if params.verbose
    T=text_waitbar(sprintf('Getting spikes and STRFs for %d cells...',length(cell_indices)));
    start_time = clock; % note when it started
end

% cycle through each cell
for cc = 1:length(cell_indices)
    
    if params.verbose
        T=text_waitbar(T,(cc)/length(cell_indices));
    end
    
    % get cell index, id
    cell_index = cell_indices(cc);
    cell_id = datarun.cell_ids(cell_index);


    % get spike times

%%%% COMMENTED OUT BY GDF ON 2014/01/09 TO FIX BUG    
    % compute spike rate at all times
%    spike_rate = histc(datarun.spikes{cell_index},datarun.triggers(1):refresh_time:movie.size*refresh_time);
    % store spikes in the relevant region
%    spikes(:,cc) = spike_rate;
%%

    % compute spike rate at all times
    spike_rate = histc(datarun.spikes{cell_index},datarun.triggers(1):refresh_time:movie.size*refresh_time+datarun.triggers(1));
    % store spikes in the relevant region
    if isempty(spike_rate)
        spikes(:,cc) = zeros(size(spikes,1),1);
    else
        spikes(:,cc) = spike_rate(1:end-1);
    end
    

    % get STA
    sta = get_sta(datarun,cell_id);
    
    % get sig stixels
    switch params.marks
        case 'marks'
            sig_stixels = datarun.stas.marks{cell_index};
        case 'simple'
            datarun = setsimplemarks(datarun, cell_id, params.simplemarkthresh);
            sig_stixels = datarun.stas.marks{cell_index};
        case 'sigstix'
            sig_stixels = significant_stixels(sta);
    end
    sig_stixels = sig_stixels(:); % Vectorize

    % identify STA frames to use
    frames = parse_frame_spec(params.frames, size(sta,4));
    nframes = length(frames);

    % Expand across color channels
    sig_stixels = repmat(sig_stixels, size(sta,3), 1);
    
    % Get desired STA frames, vectorize them
    staframes = reshape(sta(:,:,:,frames), [], nframes);
    
    % Sparsify according to sigstix
    numnonzero = sum(sig_stixels);
    strfs{cc} = spalloc(size(staframes,1), nframes, numnonzero*nframes);
    for i = 1:length(frames)
        strfs{cc}(sig_stixels,i) = staframes(sig_stixels,i);
    end
    
end


% display how long it took
if params.verbose
    fprintf('done (%0.1f seconds)\n',etime(clock,start_time));
end


% identify when to start
start_stim = floor(1+params.start_time/(movie.getRefreshTime/1000));
end_stim = floor(1+params.end_time/(movie.getRefreshTime/1000));

% when to stop
if isempty(params.stimuli)
    end_stim = movie.size;
else
    end_stim = start_stim + params.stimuli;
end


% determine if the stimulus is RGB or BW
if datarun.stimulus.independent == 't'
    isRGB = true;
else
    isRGB = false;
end



% Get generator signal values
gensigopts = {'start_stim', start_stim, 'end_stim', end_stim, 'num_frames', nframes, 'isRGB', isRGB, 'verbose', params.verbose};
gen_signals = compute_gen_signals(strfs, movie, gensigopts{:});

% pare spikes down to just the time period of interest and account for temporal duration of each RF
first_stim = start_stim + nframes - 1; % first_stim is the earliest stimulus frame for which a generator signal value can be computed
stim_frames = first_stim:end_stim;
spikes = spikes(stim_frames,:);

% compute SNLs
[fit_params, snl_params] = compute_snls(spikes, gen_signals);

% push the gensigopts onto snl_params (backward compatibility)
for i = 1:2:length(gensigopts)
    snl_params.(gensigopts{i}) = gensigopts{i+1};
end

% store in datarun
for cc = 1:length(cell_indices)
    cell_index = cell_indices(cc);
    datarun.stas.snls{cell_index}.fit_params = fit_params{cc};
    datarun.stas.snls{cell_index}.gen_signal = gen_signals(:,cc);
    datarun.stas.snls{cell_index}.spikes = spikes(:,cc);
end


