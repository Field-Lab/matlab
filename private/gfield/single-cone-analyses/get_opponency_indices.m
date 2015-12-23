function opponency_indices = get_opponency_indices(datarun, varargin)
%
% get_opponency_indices     compute opponency indices for cells in datarun
%
% usage: oppoency_indices = get_opponency_indices(datarun, varagin)
%
% arguments:           datarun - standard datarun struct 
%
% outputs:   opponency_indices - vector of opponency indices from cells in datarun
%
% optional parameters, their default values, and what they specify:
% 
% labels            []          an alternative list of cone labels -- can
%                               be used to test the impact of alternative cone mosaics, i.e. test
%                               clumping
% cell_specification  'all'     opponency indicies will only be calculated
%                               for these cells
% 
% GDF 2009-02
%


% parse inputs
p = inputParser;
p.addRequired('datarun')
p.addParamValue('cell_specification', 'all')
p.addParamValue('labels', [], @ischar)

p.parse(datarun, varargin{:})
cell_spec = p.Results.cell_specification;

% body of function
cell_indices = get_cell_indices(datarun, cell_spec);
total_cell_number = length(get_cell_indices(datarun, 'all'));

opponency_indices = zereos(total_cell_number,1);
opponency_indices(:) = NaN;  % make them NaNs instead of zeros to avoid confusion

num_RGCs = length(cell_indices);
for RGC = 1:num_RGCs
    temp_weights = datarun.cones.weights(:,cell_indices(RGC));
    temp_index = compute_opponency_index(temp_weights, datarun.cones.types);
    opponency_indices(cell_indices(RGC)) = temp_index;
end
    
