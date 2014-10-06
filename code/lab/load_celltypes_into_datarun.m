function datarun = load_celltypes_into_datarun(datarun, cell_types, verbose)
% LOAD_CELLTYPES_INTO_DATARUN
% usage: datarun = load_celltypes_into_datarun(datarun, cell_types)
%
% Abstracted from LOAD_PARAMS, LOAD_VISION_CELL_TYPES
%
% See also: LOAD_TXT_CELL_TYPES, LOAD_RRS_CELL_TYPES, ORDER_CELL_TYPES,
% STRUCT/LOAD_TXT_CELL_TYPES
%
% 2012-07 phli
%

% Save cell_types in datarun
datarun.vision.cell_types = cell_types;

% Load these cell types into datarun.cell_types if none have been defined so far
if ~isfield(datarun,'cell_types') || isempty(datarun.cell_types)
    datarun.cell_types = datarun.vision.cell_types;
    used_cell_types = true;
else
    used_cell_types = false;
end

% display which cell types were loaded
if verbose
    fprintf('\nLoaded %d cell types:\n',length(datarun.cell_types))
    show_cell_types(datarun.cell_types)
    fprintf('\n')
    
    if used_cell_types
        fprintf('These cell types are loaded into the working list of cell types\n\n')
    else
        fprintf('These cell types were NOT loaded into the working list of cell types\n\n')
    end
end
