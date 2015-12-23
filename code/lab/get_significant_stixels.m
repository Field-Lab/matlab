function returned_datarun = get_significant_stixels(datarun, cell_specification, varargin)
% get_significant_stixels     get the significant stixels from stas
%
% usage:  returned_datarun = get_significant_stixels(datarun, cell_specification, varargin)
%
% arguments:  
%             datarun - datarun structure
%  cell_specification - cells to operate on
%
% outputs:
%    returned_datarun - datarun structure
%
% optional fields in varargin, their default values, and what they specify:
%
%    verbose        false                 logical for verbose output   
%    thresh_params  []                    see significant_stixels.m
%
% gdf 2008-10-05
%   notes: gdf would like to add the ability to specify a location to find
%   stas with in datarun and the ability to specify an output field in
%   datarun
%


% SET UP OPTIONAL ARGUMENTS

p = inputParser;
p.addRequired('datarun', @isstruct);
p.addRequired('cell_specification');

p.addParamValue('thresh_params', struct, @isstruct);
p.addParamValue('verbose', false, @islogical);

% parse inputs
p.parse(datarun, cell_specification, varargin{:});
thresh_params = p.Results.thresh_params;

% verbose for clocking function
if p.Results.verbose
    fprintf('\nComputing something important...');
    start_time = clock; % note when it started
end


% logic for getting the number of cells to analyze and their indices in
% datarun
cell_indices = get_cell_indices(datarun, cell_specification);
num_cells = length(cell_indices);

    
% BODY OF THE FUNCTION

temp_datarun = datarun;
temp_datarun.stas.significant_stixels = cell(size(datarun.stas.stas));
for cll = 1:num_cells
    temp_sig_stix = significant_stixels(datarun.stas.stas{cell_indices(cll)}, thresh_params);
    temp_datarun.stas.significant_stixels{cell_indices(cll)} = temp_sig_stix;
end

returned_datarun = temp_datarun;

% display how long it took
if p.Results.verbose
    fprintf(' done (%0.1f seconds)\n',etime(clock,start_time));
end