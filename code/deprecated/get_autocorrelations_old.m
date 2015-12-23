function datarun = get_autocorrelations(datarun, cell_spec, varargin)
% MY_FUNCTION    compute the autocorrelation for all cells specified in
%                cell_spec
%
% usage:  datarun = get_auto_correlations(datarun, cell_spec, varargin)
%
% arguments:  datarun - datarun struct 
%           cell_spec - which cells (see get_cell_indices for options)
%            varargin - struct or list of optional parameters (see below)
%
% outputs:    datarun - datarun struct with results stored in field
%                       autocorrelation.probabilities{}
%                       autocorrelation.bins{}    
%
% optional parameters, their default values, and what they specify:
%
%
% bin_size       0.001              bin size for autocorrelogram (in
%                                   seconds)
% duration        0.1               duration to calculate autocorrelation
%                                   (in seconds)       
%
%
%
% 2011-02  GDF
%


% SET UP OPTIONAL ARGUMENTS

p = inputParser;

% specify list of optional parameters
p.addParamValue('bin_size', 0.001, @isnumeric);
p.addParamValue('duration', 0.1, @isnumeric);

% resolve user input and default values
p.parse(varargin{:});

% BODY OF FUNCTION

% get number of cells and indices
cell_indices = get_cell_indices(datarun, cell_spec);
num_cells = length(cell_indices);

% initialize autocorrelation field if it does not already exist
if ~isfield(datarun, 'autocorrelation')
    datarun.autocorrelation = cell(length(datarun.cell_ids),1);
end

% get autocorrelation and bins for cells of interest
for cc = 1:num_cells
    spike_times = datarun.spikes{cell_indices(cc)};
    [probabilities, bins] = autocorrelation(spike_times, p.Results.bin_size, p.Results.duration);

    % store information in datarun
    datarun.autocorrelation{cell_indices(cc)}.probabilities = probabilities;
    datarun.autocorrelation{cell_indices(cc)}.bins = bins;
end

    