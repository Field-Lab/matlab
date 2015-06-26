function [type_name,type_index]=find_cell_type(datarun, cellID)

% cellID for vision cell ID, not datarun cell index

for i=1:length(datarun.cell_types)
    if ~isempty(datarun.cell_types{i}.cell_ids)
        tmp=find(datarun.cell_types{i}.cell_ids==cellID, 1);
        if ~isempty(tmp)
            type_name = datarun.cell_types{i}.name;
            type_index = i;
            break;
        end
    end
end
    