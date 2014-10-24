function type = getCellType(cell_types,rgc_id)

% input: cell array of cone types, index of an rgc
% output: string containing the type of rgc

rgc_id = rgc_id(:);
nTypes = length(cell_types);
nCells = size(rgc_id,1);
type = cell(nCells,1);
for iCell=1:nCells
    found = 0;
    for iType=1:nTypes
        if sum(cell_types{iType}.cell_ids == rgc_id(iCell));
            type{iCell} = cell_types{iType}.name;
            found = 1;
        end
    end
    if ~found
        fprintf('cannont find type for cell %g\r',rgc_id(iCell));
    end
end
