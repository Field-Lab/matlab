function returned_datarun = get_stas_color_transformed(datarun, cell_specification, transform, varargin)
% get_stas_color_transformed     Applies transform to color dimension of stas
%
% usage:  datarun = get_stas_color_transformed(sta, transform, params)
%
% arguments:  datarun - datarun struct 
%           transform - a 3x3 matrix 
%
% outputs:     
%    returned_datarun - color transformed stas are returned in
%                       returned_datarun.stas.color_transformed_stas{}
% optional fields in varargin, their default values, and what they specify:
%
%    cell_ids       []                    numerical list of cell ids
%    verbose        false                     logical for verbose output   
%
%
% GDF 2008-10-2
%

% SET UP OPTIONAL ARGUMENTS

p = inputParser;
p.addRequired('datarun', @isstruct);
p.addRequired('cell_specification')
p.addRequired('transform', @isnumeric);
p.addParamValue('verbose', false, @islogical)

% parse inputs
p.parse(datarun, cell_specification, transform, varargin{:});
verbose = p.Results.verbose;

if verbose
    fprintf('\nComputing something important...');
    start_time = clock; % note when it started
end


% get indexes for cell IDs
cell_indices = get_cell_indices(datarun, cell_specification);
cell_number = length(cell_indices);

% APPLY TRANSFORM TO STA 

temp_datarun = datarun;
temp_datarun.stas.color_transformed_stas = cell(size(datarun.stas.stas));
for cll = 1:cell_number
    temp_sta = sta_color_transformed(datarun.stas.stas{cell_indices(cll)}, transform, 'verbose', verbose);
    temp_datarun.stas.color_transformed_stas{cell_indices(cll)} = temp_sta;
end

returned_datarun = temp_datarun;

% display how long it took
if verbose
    fprintf(' done (%0.1f seconds)\n',etime(clock,start_time));
end