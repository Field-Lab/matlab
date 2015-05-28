function [datarun] = load_data(varargin)
% LOAD_DATA         Initialize datarun variable
%
% usage:  datarun = load_data( <spec>, params )
%         datarun = load_data( { {<spec>}, {<spec>}, ..., {<spec>} } params  )
%
% <spec> can specify a datarun in one of three ways:
%
%  <spec> = experiment, condition
%           experiment is a string specifying the piece, e.g. '2006-07-05-1'
%           condition is a string specifying the series, e.g. 'rf-0-mg'
%
%  <spec> = rrs_prefix
%           string specifying the path for vision files,
%               e.g. '/Analysis/Greschner/2006-07-05-1/data000/data000'
%           if the first character is not '/', the path to the server 
%               (global server_path) is prepended
%
%  <spec> = datarun
%           struct of an existing datarun
%           based on the entries in .names the acording data are loaded and
%           added, while the rest of the struct is unchanged
%
% When one specification is provided, load_data returns one datarun struct.
% When multiple specifications are provided, load_data returns a cell array
% of datarun structs.
%
%
% optional fields in params, their default values, and what they specify:
%
% verbose               false   prints datarun information as it is loaded
% load_index            true	load index file information
%                                   NOTE: if <spec> = rrs_prefix, load_index defaults to false
% load_neurons          false   load neurons file
% load_ei               false 	load ei file
% load_params           false 	load params file
% load_sta           	false  	load sta file
% load_obvius_sta_fits	false   load obvius sta fits
% load_*_params         []      struct of arguments for load_*, where * = 
%                                   {index, neurons, ei, params, sta, obvius_sta_fits}
% load_all              false   sets load_* to true, unless otherwise specified
% set_polarities        {}      if not empty, passes as args to set_polarities
% map                   false   map cell_types for 2 dataruns
%
%
% usage examples:
%
% datarun = load_data(datarun);
%
% datarun = load_data('2007-07-05-1','rf-0-mg');
%
% datarun = load_data('/Analysis/Greschner/2007-07-05-1/data000/data000');
%
% datarun = load_data({{'2007-07-05-1','rf-0-mg'},...
%               {'2007-07-05-1','rf-1-mapped-0-mg'}});
%
% datarun = load_data(datarun, struct('verbose',true));
%
%
% gauthier 2008-03
% greschner change verbose output at % LOAD ALL PARTS, AS POSSIBLE/SPECIFIED
% greschner add map cell_types 
% greschner fix 'obvius_fit_path'


%   IMPLEMENT THESE OPTIONS
%
%    datarun{1}.names.rrs_params_path='2007-02-06-3/data009-gdf/data009.params';
%    datarun{1}.stimulus.stixel_width=16;
%    datarun{1}.stimulus.stixel_height=datarun{1}.stimulus.stixel_width;
%    datarun{2}.names.rrs_neurons_path='2007-02-06-3/data013-tb/data013-from-data009/data013-map.neurons';
%    datarun{2}.stimulus.type='gray';
%    opt=struct('verbose',true,'map',true);
%    datarun=load_data(datarun,opt);




% determine how many dataruns were specified, and load them one at at time
% with the internal function load_one_datarun


% if dataruns were specified in a cell array, there could be several, so load each in turn

% if the first argument is a cell array ...
if iscell(varargin{1})
    % there should be only one cell array, potentially followed by a params struct.
    % so return an error if there are more than 2 total arguments, or exactly 2 and the second argument isn't the params struct
    if length(varargin) > 2 || (length(varargin)==2 && ~isstruct(varargin{2}))
        error('Type HELP LOAD_DATA for appropriate usages.')
    else
        % otherwise, assume the cell array is a list of datarun specifications
        
        % get list of dataruns
        datarun_list = varargin{1};
        
        % get params struct from input arguments, if it exists
        if length(varargin) == 2
            params = varargin{2};
        else % otherwise, set up an empty struct
            params = struct;
        end

        % load each datarun in turn
        for dd = 1:length(datarun_list)
            this_datarun = datarun_list{dd};
            
            % if it's a cell array of strings...
            if iscellstr(this_datarun)
                % pass them in separately
                datarun{dd} = load_one_datarun(this_datarun{:},params);
                
                % otherwise, assume it's a datarun struct
            elseif isstruct(this_datarun)
                datarun{dd} = load_one_datarun(this_datarun,params);
            else
                
                % if not give an error
                error('Type HELP LOAD_DATA for appropriate usages.')
            end
        end
        
        % if map is true, load "mapped cell types" into the working list of cell types for the first two dataruns
        if length(varargin) > 1 && isstruct(varargin{end}) && isfield(varargin{end},'map') && varargin{end}.map
            datarun{1}.cell_types=map_cell_types(datarun{1}, datarun{2}, struct('verbose',true));
            datarun{2}.cell_types=datarun{1}.cell_types;
        end

    end

else  % otherwise, assume there is only one datarun
    datarun = load_one_datarun(varargin{:});
end


function pth = infer_file_spec(pth)
splitpath = strsplit(pth, '/');

% If no filespec was given, assume it matches the directory
if length(splitpath) == 2
    splitpath{end+1} = splitpath{end};
    pth = strjoin(splitpath, '/');
end


function datarun = load_one_datarun(varargin)
% datarun specification arguments must be match one of these 6 forms:
%
% <struct>                          datarun struct
% <string>                          rrs prefix
% <string>, <string>                experiment & condition
%  
% <struct>, <struct>                datarun struct, plus params struct
% <string>, <struct>                rrs prefix, plus params struct
% <string>, <string>, <struct>      experiment & condition, plus params struct



% PEEL OFF THE PARAMS STRUCT, IF IT EXISTS

% check to be sure there aren't too many arguments
if length(varargin) > 3 || (length(varargin) > 2 && ~isstruct(varargin{end}))
    error('Type HELP LOAD_DATA for appropriate usages.')
    
else  % if there is params struct, separate it from the other arguments
    if length(varargin) > 1 && isstruct(varargin{end})
        % put the params struct into one variable
        params = varargin{end};
        % and the datarun specification argument(s) into another
        datarun_spec = {varargin{1:end-1}};
    else
        % if there is no params struct, then all arguments are the datarun specification 
        datarun_spec = varargin;
        params = struct;
    end
end


% PREPARE datarun STRUCT FOR THE LOADING FUNCTIONS
% if it's not a datarun struct, then create a datarun struct with the appropriate fields

% determine which kind of specification was used
switch length(datarun_spec)
    
    case 1 % a struct or vision file path
        switch class(datarun_spec{1})
            
            
            case 'struct' % datarun struct
                datarun = datarun_spec{1};
                
                
            case 'char' % file path
                % this could be an rrs prefix, or a single file (the latter is detected by the presence of a period)
                
                input_path = datarun_spec{1};
                
                % If there is only a directory and no filespec, assume
                % filespec matches directory.  I.e.,
                % "2010-07-22-0/data000" ==> "2010-07-22-0/data000/data000"
                input_path = infer_file_spec(input_path);

                % if the file path does not start with a '/', add server path to the beginning
                input_path = add_base_path(input_path, server_path());
                
                % if the path ends in a file extension, then only load that particular file.
                % otherwise, if there is no extension, treat it as an rrs prefix
                
                % look for a period
                period_index = strfind(input_path,'.');                
                if isempty(period_index)
                    % if there is no period, assume this is an rrs prefix
                    datarun.names.rrs_prefix = input_path;
                else
                    % otherwise assume there's a file extension after the period

                    % only look at the final period (in case there are several)
                    period_index = period_index(end);

                    % be sure this period wasn't the last character
                    if period_index == length(input_path)
                        error('File path specification should not end in a period.')
                    end

                    % look at the file extension...
                    file_extension = input_path(period_index+1:end);
                    
                    switch file_extension
                        case {'params','neurons','ei','sta','globals'}
                            % ...and if it is recognized, load that file

                            % if the file exists...
                            if exist(input_path,'file')
                                % ...set flag and field to load it
                                params.(['load_' file_extension ])= true;
                                datarun.names.(['rrs_' file_extension '_path']) = input_path;
                                
                                % also set up the prefix
                                datarun.names.rrs_prefix = input_path(1:period_index-1);
                            else
                                % otherwise give an error
                                error('No file found at ''%s'' ',input_path)
                            end

                        otherwise
                            % if the extension was not recognized, give an error
                            error('File extension ''%s'' is not a valid file for loading.',file_extension)
                    end
                end

            otherwise % error
                error('Type HELP LOAD_DATA for appropriate usages.')
        end
        
    case 2 % experiment and condition
        if ischar(datarun_spec{1}) && ischar(datarun_spec{2})
            datarun.names.experiment = datarun_spec{1};
            datarun.names.condition = datarun_spec{2};
        else
            error('Type HELP LOAD_DATA for appropriate usages.')
        end
        
    otherwise % error
        error('Type HELP LOAD_DATA for appropriate usages.')
end





% SET UP OPTIONAL PARAMETERS


% if params.load_all is true, set all 'load_*' fields in params to true,
% but don't overwrite any fields which already exist
if isfield(params,'load_all') && params.load_all
    % set list of options to mark as true
    load_options = {'load_index','load_neurons','load_ei','load_params','load_sta','load_obvius_sta_fits'};
    % go through list of load options
    for oo = 1:length(load_options)
        % if there is no specified value
        if ~isfield(params,load_options{oo})
            % set it to true
            params.(load_options{oo}) = true;
        end
    end
end


% specify default parameters

defaults.verbose = false;
defaults.load_all = false;

% default to loading the index file only if the appropriate fields exist
if isfield(datarun.names,'experiment') && isfield(datarun.names,'condition')
    defaults.load_index = true;
else defaults.load_index = false; end

defaults.load_index_params = struct;

defaults.load_neurons = false;
defaults.load_neurons_params = struct;

defaults.load_ei = false;
defaults.load_ei_params = struct;

defaults.load_params = false;
defaults.load_params_params = struct;

defaults.load_sta = false;
defaults.load_sta_params = struct;

defaults.load_obvius_sta_fits = false;
defaults.load_obvius_sta_fits_params = struct;

defaults.set_polarities = {};
defaults.map = false;

% combine user and default parameters
params = default_params(defaults, params);




% LOAD ALL PARTS, AS POSSIBLE/SPECIFIED

% load index file
if params.load_index
    % only load if required fields exist
    if isfield(datarun.names,'experiment') && isfield(datarun.names,'condition')
        if params.verbose; disp('load index file'); end;
        datarun = load_index(datarun, params.load_index_params);
    else
        if params.verbose; disp('skipping index file'); end;
    end
end


% expand rrs_prefix - creates datarun.names.rrs_neurons_path, etc. based on
% names.rrs_prefix
datarun = expand_rrs_prefix(datarun);

% fill in obvius fit path
if ~isfield(datarun.names,'obvius_fit_path')
    if isfield(datarun.names,'experiment') && isfield(datarun.names,'condition')
        datarun.names.obvius_fit_path =  [server_path datarun.names.experiment '/' datarun.names.condition '/'];
    else
        datarun.names.obvius_fit_path = [];
    end
end  

% assign short name
datarun.names.short_name = short_name_from_names(datarun.names);

% load neurons file
if params.load_neurons    
    rrs_neurons_path = get_rrs_neurons_path(datarun);

    % only load if required field exists
    if ~isempty(rrs_neurons_path)
        datarun.names.rrs_neurons_path = add_server_path(rrs_neurons_path);
        if params.verbose
            disp(sprintf('load neurons file: %s',datarun.names.rrs_neurons_path));
        end
        datarun = load_neurons(datarun, params.load_neurons_params);
    else
        if params.verbose; disp('skipping neurons file'); end;
    end
end


% load ei file
if params.load_ei
    % only load if required field exists
    if isfield(datarun.names,'rrs_ei_path') & isfield(datarun,'cell_ids')
        datarun.names.rrs_ei_path = add_server_path(datarun.names.rrs_ei_path);
        if params.verbose;
            disp(sprintf('load ei file: %s',datarun.names.rrs_ei_path));
        end;
        datarun = load_ei(datarun, datarun.cell_ids, params.load_ei_params);
    else
        if params.verbose; disp('skipping ei file - neurons file is required'); end;
    end
end


% load params file
if params.load_params
    % only load if required field exists
    if isfield(datarun.names,'rrs_params_path')
        datarun.names.rrs_params_path = add_server_path(datarun.names.rrs_params_path);
        if params.verbose;
            disp(sprintf('load params file: %s',datarun.names.rrs_params_path));
        end;
        datarun = load_params(datarun, params.load_params_params);
    else
        if params.verbose; disp('skipping params file'); end;
    end
end


% load sta file
if params.load_sta
    % only load if required field exists
    if isfield(datarun.names,'rrs_sta_path')
        datarun.names.rrs_sta_path = add_server_path(datarun.names.rrs_sta_path);
        if params.verbose;
            disp('load sta file');
            disp(sprintf('load sta file: %s',datarun.names.rrs_sta_path));
        end;
        datarun = load_sta(datarun, params.load_sta_params);
    elseif isfield(datarun.names, 'split_rrs_sta_path')
        if params.verbose;
            disp('load sta file');
            disp(sprintf('load sta file: %s',datarun.names.split_rrs_sta_path));
        end;
        datarun = load_sta(datarun, params.load_sta_params);
        if params.verbose; disp('skipping sta file'); end;
    end
end


% load obvius sta fits
if params.load_obvius_sta_fits
    if params.verbose; 
        disp('load obvius sta fits'); 
    end;
    %if isfield(datarun.names,'obvius_fit_path')
    if ~isempty(datarun.names.obvius_fit_path)
        %datarun = load_obvius_sta_fits(datarun, params.load_obvius_sta_fits_params);
        datarun = load_obvius_sta_fits(datarun);
    else
        if params.verbose; disp('skipping obvius sta fits'); end;
    end
end


% Set polarities?
if ~isempty(params.set_polarities)
    datarun = set_polarities(datarun, params.set_polarities{:});
end