function n = getCellTypeNum(cell_types,cell_type)

% input: cell array of cell types, desired cell type
% output: number of cells of the desired type

nTypes = length(cell_types);

found = 0;
for iType=1:nTypes
    if strcmpi(cell_types{iType}.name,cell_type)
        n = length(cell_types{iType}.cell_ids);
        found = 1;
        break
    end
end

if found == 0;
    error('(getCellTypeNum) cell type %s not found\r',cell_type);
end