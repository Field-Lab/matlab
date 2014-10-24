function ids = getTypeIds(cell_types,thisType)

% input: array of cell types
% output: ids of all cells having that type

nTypes = length(cell_types);

found = 0;
for iType=1:nTypes
    if strcmpi(cell_types{iType}.name,thisType)
        ids = cell_types{iType}.cell_ids;
        found = 1;
    end
end
if ~found
    fprintf('(getTypeIds) cannot find type\r');
    ids = [];
end
