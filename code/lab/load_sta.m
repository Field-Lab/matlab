function datarun = load_sta(datarun, varargin)
% LOAD_STA     Load information from a sta file
%
% usage:  datarun = load_sta(datarun, params)
%
% arguments:  datarun - datarun struct with these fields
%                       	datarun.names.rrs_sta_path - absolute path to sta file,
%                               e.g. '/Analysis/Greschner/2005-04-26-1/data000/data000.sta'
%                           datarun.cell_ids - list of cell IDs in the datarun
%              params - struct of optional parameters (see below)
%
% outputs:    datarun - datarun struct with the following field added or updated
%
%     	datarun.stas.stas       cell array of 4-d matrices
%                                   dimensions are height, width, color, and frames
%                                   elements are single precision floats
%
%
% optional fields in params, their default values, and what they specify:
%
% load_sta            	'all'	which cells to load an STA for (for other acceptable values,
%                                   see get_cell_indices and examples below)
% verbose            	false	print information about loading
% show_mem              false   after loading each cell, show how much memory datarun occupies
% frames                ':'     which frames to load
%                                   ':' - all
%                                     n - frame (-n+1) through frame 0
% sync_stimulus         true	if datarun.stimulus doesn't match the STAs, give an error
%                               if datarun.stimulus doesn't exist, set it to match the STAs
% sync_cell_ids         true    same as sync_stimulus, applied to list of cell IDs
% save_sta              true    put STAs into the datarun struct
% keep_java_sta         true    put the java sta object into datarun.stas.java_sta
% save_rf               false   compute and save summary frame (see rf_from_sta for options)
% rf_params          	[]     	params struct for rf_from_sta (only used if save_rf is true)
% guess_stimulus        true    fill in missing parameters of the stimulus
%                                   calls function 'guess_stimulus'
%
%
%
% examples:
%
%	datarun = load_sta(datarun);
%
%	datarun = load_sta(datarun, struct('load_sta',[1153 3005]));
%
%	datarun = load_sta(datarun, struct('load_sta','ON parasol'));
%
%	datarun = load_sta(datarun, struct('load_sta',{{1,2}}));
%       note that matlab syntax requires TWO sets of curly brackets to set
%       the value of a field in a struct to be a cell array
%
%   datarun = load_sta(datarun, struct('load_sta',[]));
%   	this will synchronize datarun.stimulus and the STAs without loading anything
%
%
%
% 2008-10  gauthier
% 2012-10  phli, loads STA fits if possible
%
%



% SET UP OPTIONAL ARGUMENTS

p = inputParser;
p.addRequired('datarun', @isstruct);
p.addOptional('java_sta', [],@isjava);
% specify list of optional parameters
p.addParamValue('load_sta', 'all');
p.addParamValue('verbose', false);
p.addParamValue('show_mem', false);
p.addParamValue('frames', ':');
p.addParamValue('sync_stimulus', true);
p.addParamValue('sync_cell_ids', true);
p.addParamValue('save_sta', true);
p.addParamValue('keep_java_sta', true);
p.addParamValue('save_rf', false);
p.addParamValue('rf_params', struct);
p.addParamValue('guess_stimulus', true);

% resolve user input and default values
p.parse(datarun, varargin{:});

% get params struct
params = p.Results;




% ensure java stas object is loaded

% use the user argument, if provided
if ~isempty(params.java_sta)
    java_sta = params.java_sta;
else
    % otherwise, check if it's already loaded
    if isfield(datarun,'stas') && isfield(datarun.stas,'java_sta') && ~isempty(datarun.stas.java_sta)
        java_sta = datarun.stas.java_sta;
    else
        % if not, load it
        java_sta = load_java_sta(datarun,'verbose',params.verbose);
    end
end

datarun.stas.depth = java_sta.getSTADepth();

% get stimulus parameters from the STAs file
stimulus_from_stas = stimulus_from_java_stas(java_sta);

% get additional params from globals file
stim_from_globals = struct();
if isfield(datarun.names,'rrs_globals_path')
    datarun.names.rrs_globals_path = add_server_path(datarun.names.rrs_globals_path);
end
datarun = load_globals(datarun);
if isfield(datarun, 'globals') && ~isempty(datarun.globals)
    stim_from_globals = stimulus_from_globals(datarun.globals);
end

if params.sync_stimulus
    % if datarun.stimulus doesn't match the STAs, give an error.
    % if datarun.stimulus doesn't exist, define it to match the STAs
    datarun = sync_stimulus(datarun, stimulus_from_stas, 'STAs');
    
    % globals overrides STAs for interval, since the STAs logic doesn't
    % know about OLED or other non-120-Hz stim sources
    if isfield(stim_from_globals, 'interval') && stim_from_globals.interval > 0
        datarun.stimulus.interval = stim_from_globals.interval;
    end
    
    % Remaining params should sync between STAs and globals or should be
    % new ones that we can safely bring in from globals.
    datarun = sync_stimulus(datarun, stim_from_globals, 'Globals file');
    
    fprintf('\nStimulus information synchronized.\n')
end

% fill in missing stimulus parameters, if desired
if params.guess_stimulus
    datarun.stimulus = guess_stimulus(datarun.stimulus);
end



% get cell IDs from STAs file
cell_ids_from_stas = sort(java_sta.getIDList)';

% if datarun.cell_ids doesn't match the STAs, give an error.
% if datarun.cell_ids doesn't exist, define it to match the STAs
if params.sync_cell_ids
    if isfield(datarun.names, 'rrs_sta_path')
        datarun = sync_cell_ids(datarun, cell_ids_from_stas,  sprintf('sta file %s',datarun.names.rrs_sta_path));
    elseif isfield(datarun.names, 'split_rrs_sta_path')
        datarun = sync_cell_ids(datarun, cell_ids_from_stas,  sprintf('sta file %s',datarun.names.split_rrs_sta_path));;
    end
        
    fprintf('\nCell ids synchronized.\n')
end




% get total number of cells
cell_count = length(cell_ids_from_stas);

% make fields for stas, rfs, etc in advance (even if they won't be loaded)
datarun = make_fields(datarun,cell_count,{'stas','rfs','marks','rf_coms','time_courses'});



% if STAs or RFs will be saved...
if params.save_sta || params.save_rf


    % get list of cell IDs to load
    cell_nums = get_cell_indices(datarun, params.load_sta);


    % display how many cells will be loaded
    if params.verbose;
        fprintf('\nLoading STAs for %d cells ',length(cell_nums));
        start_loading = clock; % note when it started
    end


    % go through list of cells and load up each STA
    for cell_num = cell_nums

        % show a tick of progress
        if params.verbose;fprintf('.');end

        % get STA from the java object
        sta = sta_from_java_stas(java_sta, datarun.cell_ids(cell_num), params.frames, stimulus_from_stas.independent);

        % save STA, if desired
        if params.save_sta
            datarun.stas.stas{cell_num} = sta;
        end

        % compute RF, if desired
        if params.save_rf
            datarun.stas.rfs{cell_num} = rf_from_sta(sta, params.rf_params);
        end


        % show datarun size
        if params.show_mem
            a=whos;
            b=whos('datarun');
            fprintf('\n%0.1f MB in datarun, %0.1f MB otherwise (%0.1f MB total)',b.bytes/(2^20), sum([a.bytes])/(2^20)-b.bytes/(2^20), sum([a.bytes])/(2^20))
        end
    end


    % display how long it took
    if params.verbose;
        fprintf('\n     done (%0.1f minutes)\n',etime(clock,start_loading)/60);
    end

end

% If sta fits are available from only one place, go ahead and use those
if isfield(datarun, 'vision') && isfield(datarun.vision, 'sta_fits') && (~isfield(datarun, 'obvius') || ~isfield(datarun.obvius, 'sta_fits'))
    datarun = get_sta_fits_from_vision(datarun);
elseif isfield(datarun, 'obvius') && isfield(datarun.obvius, 'sta_fits') && (~isfield(datarun, 'vision') || ~isfield(datarun.vision, 'sta_fits'))
    datarun = get_sta_fits_from_obvius(datarun);
end

% keep the java sta object, or close it
if params.keep_java_sta
    datarun.stas.java_sta = java_sta;
else
    java_sta.close
end



function datarun = make_fields(datarun,cell_count,field_list)

% ensure datarun.stas exists
if ~isfield(datarun,'stas')
    datarun.stas = struct;
end

% for each field
for ff = 1:length(field_list)
    % ensure the field exists...
    if ~isfield(datarun.stas,field_list{ff})
        datarun.stas.(field_list{ff}) = cell(cell_count,1);
    else
        % ... and is the correct length
        if length(datarun.stas.(field_list{ff})) < cell_count
            datarun.stas.(field_list{ff}){cell_count} = [];
        end
    end
end

