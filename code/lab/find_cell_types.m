function cell_type_numbers = find_cell_types(datarun, cell_spec,params)
% FIND_CELL_TYPES     identify the cell types of a list of cells
%
% usage:  cell_type_numbers = find_cell_types(datarun, cell_spec, params)
%
% arguments:  datarun - datarun struct
%           cell_spec - see get_cell_indices for options
%              params - struct of optional arguments, see below
%
% outputs:   cell_type_numbers - NxM vector of cell type numbers
%
%
% optional fields in params, their default values, and what they specify:
%
% return            'first'  	if a cell occurs in multiple types, which ones to return
%                                   'first' - only return the number of the first listed cell type
%                                               in which the cell appears.  With this option,
%                                               the result is a vector.
%                                   'all'   - return the numbers of all cell types
%                                               in which the cell appears.  With this option,
%                                               the result is a NxM matrix where each row is
%                                               list of cell types in which the cell appears.
%                                               For cells appearing in fewer than the maximum 
%                                               number of types, the row is padded with zeros.
%
%
% 2008-09  gauthier
%
%


% SET UP OPTIONAL ARGUMENTS

% if not specified, make params empty
if ~exist('params','var');params = [];end

% specify default parameters
defaults.return = 'first';

% combine user and default parameters
params = default_params( defaults, params);



% BODY OF THE FUNCTION

% get cell indices
cell_indices = get_cell_indices(datarun,cell_spec);


% get total number of cell types
if isfield(datarun,'cell_types')
    num_cell_types = length(datarun.cell_types);
else
    % if cell types are not defined, return empty
    cell_type_numbers = [];
    return
end
    

% store lists of cell ids in a matrix.  each row lists the cell ids in a single cell type.
% matrix is padded with zeros.  if a cell type has no cell ids, the row is all zeros.

% initialize matrix
cell_id_lists = [];

% go through each cell type
for cc = 1:num_cell_types

    if isempty(datarun.cell_types{cc})
        % if a cell type is empty, fill in the row with zeros
        cell_id_lists(cc,1) = 0;
    else
        % otherwise get the list of cell ids
        list = datarun.cell_types{cc}.cell_ids;
        cell_id_lists(cc,1:length(list)) = list;
    end
end


% convert list of cell numbers to cell ids
cell_ids = datarun.cell_ids(cell_indices);



% go through list of cells
for cc = 1:length(cell_ids)
    
    % note all the cell types in which the cell was found
    cell_found = find(max(cell_id_lists == cell_ids(cc),[],2));
    
    if isempty(cell_found)
        % if the cell was not found, return 0
        cell_type_numbers(cc,1) = 0;
    else
        % if the cell was found, return either...
        switch params.return
            case 'first'
                % the first cell type in which it appeared
                cell_type_numbers(cc,1) = cell_found(1);
            case 'all'
                % or all cell types in which it appeared
                cell_type_numbers(cc,1:length(cell_found)) = cell_found;
        end
    end   
end
    
