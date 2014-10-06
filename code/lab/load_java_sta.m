function java_sta = load_java_sta(datarun, varargin)
% load_java_sta     Ensure the java_stas object is in datarun
%
% usage:  datarun = load_java_sta(datarun, varargin)
%
% arguments:     datarun - datarun struct
%
% outputs:       datarun - datarun struct with field datarun.stas.java_sta
%
%
% optional params, their default values, and what they specify:
%
% verbose           false               show output
%
%
% gauthier 2008-10
%



% SET UP OPTIONAL ARGUMENTS

p = inputParser;
p.addRequired('datarun',@isstruct)
% specify list of optional parameters
p.addParamValue('verbose', false);

% resolve user input and default values
p.parse(datarun, varargin{:});

% get params struct
params = p.Results;




% BODY OF THE FUNCTION

% if java_stas is not already loaded
if ~isfield(datarun,'stas') || ~isfield(datarun.stas,'java_stas')
    
    % show what's going on
    if params.verbose;
        fprintf('\nLoading Java STA object...');
        start_loading = clock; % note when it started
    end
    
    % If STAFile has been created recently, after we ran the original
    % LOAD_DATA for datarun, need to rerun EXPAND_RRS_PREFIX
    if ~isfield(datarun.names, 'rrs_sta_path') && ~isfield(datarun.names, 'split_rrs_sta_path')
        datarun = expand_rrs_prefix(datarun, 'quiet', true);
    end
    
    % load object
    if isfield(datarun.names, 'rrs_sta_path')
        java_sta = edu.ucsc.neurobiology.vision.io.STAFile(datarun.names.rrs_sta_path);
    elseif isfield(datarun.names, 'split_rrs_sta_path')
        java_sta = edu.ucsc.neurobiology.vision.io.SplitSTAFile.fromBaseFilename(datarun.names.split_rrs_sta_path);
    else
        error('No rrs_sta_path or split_rrs_sta_path in datarun.names for %s', datarun.names.short_name);
    end

    % display how long it took
    if params.verbose;
        fprintf(' done (%0.1f seconds)\n',etime(clock,start_loading));
    end
end
