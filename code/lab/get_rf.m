function rf = get_rf(datarun,cell_id,varargin)
% get_rf     return the RF for a cell.  If the STA is not loaded, call get_sta
%
% usage:  rf = get_rf(datarun,cell_id,varargin)
%
% arguments:  datarun - datarun struct
%             cell_id - cell id
%            varargin - struct or list of optional parameters (see below)
%
% outputs:     rf - YxXxC matrix of the RF
%
%
% optional params, their default values, and what they specify:
%
% where     'rfs'       name of field in datarun.stas where the rfs are located
% polarity  false       multiply the returned value by the polarity from datarun.stas.polarities
%
%
% gauthier 2008-10
%


% SET UP OPTIONAL ARGUMENTS

p = inputParser;

% specify list of optional parameters
p.addParamValue('where', 'rfs',@ischar);
p.addParamValue('polarity', false,@islogical);
p.addParamValue('color_transform', [], @isnumeric)

% resolve user input and default values
p.parse(varargin{:});

% get params struct
params = p.Results;






% BODY OF THE FUNCTION

% get cell index
cell_index = get_cell_indices(datarun,cell_id);


% if the rf is present in datarun, grab it
if isfield(datarun,'stas') && isfield(datarun.stas,params.where) && ~isempty(datarun.stas.(params.where){cell_index})
    rf = datarun.stas.(params.where){cell_index};
else
    % if the RF is not present, get the STA, and then the RF
    sta = get_sta(datarun,cell_id);
    % use marks, if present
    if isfield(datarun,'stas') && isfield(datarun.stas,'marks') && ~isempty(datarun.stas.marks{cell_index})
        rf = rf_from_sta(sta,'sig_stixels',datarun.stas.marks{cell_index});
    else
        rf = rf_from_sta(sta);
    end
end


% multiply by polarity if desired
if params.polarity && ...
        isfield(datarun,'stas') && isfield(datarun.stas,'polarities') && ...
        length(datarun.stas.polarities) >= cell_index && ...
        ~isempty(datarun.stas.polarities{cell_index}) && ...
        isscalar(datarun.stas.polarities{cell_index}) && ...
        datarun.stas.polarities{cell_index} ~= 0
    rf = rf * datarun.stas.polarities{cell_index};
end


% apply color transform if specified
if ~isempty(params.color_transform);
    [num_rows, num_cols, num_pages] = size(rf);
    reshaped_rf = reshape(rf, [], num_pages);
    transformed_rf = reshaped_rf * params.color_transform;
    rf = reshape(transformed_rf, num_rows, num_cols, []);
end