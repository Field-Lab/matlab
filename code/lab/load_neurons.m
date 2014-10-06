function datarun = load_neurons(datarun, varargin)
% LOAD_NEURONS     Load information from a neurons file
%
% usage:  datarun = load_neurons(datarun, <params>)
%
% arguments:  datarun - datarun struct with field specifying the neurons file path
%                         (datarun.names.rrs_neurons_path - absolute path to neurons file,
%                          e.g. '/Analysis/Greschner/2005-04-26-1/data000/data000.neurons')
%              params - struct or list of optional parameters (see below)
%
% outputs:    datarun - datarun struct with the following fields added, as possible
%
% 	piece information:
%     	datarun.cell_ids
%     	datarun.spikes
%     	datarun.channels
%     	datarun.triggers
%     	datarun.duration
%     	datarun.sampling_rate
%
% NOTE: datarun.sampling_rate is automatically set to 20000, because the low level function load_rrs_neurons
%       can not accept other sampling rates.
%
%
% optional fields in params, their default values, and what they specify:
%
% load_spikes    	'all'       load spike times for specified cells (see get_cell_indices for specification options)
% sync_cell_ids    	true        ensure list of cell ids matches what's already in datarun
% sort_cell_ids   	false       sort cell ids
%
%
%
% 2008     gauthier
% 2009-09  gauthier, added sort_cell_ids option
%
%



% SET UP OPTIONAL ARGUMENTS

p = inputParser;

% specify list of optional parameters
p.addParamValue('load_spikes', 'all');
p.addParamValue('sync_cell_ids', true);
p.addParamValue('sort_cell_ids', false);

% resolve user input and default values
p.parse(varargin{:});

% get params struct
params = p.Results;



% CALL A LOW LEVEL FUNCTION
[spikes, extras] = load_rrs_neurons(datarun.names.rrs_neurons_path, params.load_spikes);


% STORE IN DATARUN
datarun.triggers = extras.triggers;
datarun.duration = extras.duration;
datarun.sampling_rate = 20000;

% sort cell IDs, if desired
if params.sort_cell_ids
    
    % make new variable to store sorted list
    cell_ids = zeros(size(extras.cell_ids));
    % identify sort order
    [junk, sort_order] = sort(extras.cell_ids);
    % arrange in sorted order
    for ss=1:length(sort_order)
        cell_ids(ss) = extras.cell_ids(sort_order(ss));
        datarun.channels(ss) = extras.channels(sort_order(ss));
        datarun.spikes{ss} = spikes{sort_order(ss)};
    end
else
    cell_ids = extras.cell_ids;
    datarun.channels = extras.channels;
    datarun.spikes = spikes;
end

% synchronize cell IDs, if desired
if params.sync_cell_ids
   datarun = sync_cell_ids(datarun, cell_ids', sprintf('neurons file %s',datarun.names.rrs_neurons_path));
end

datarun = build_cell_nums(datarun);