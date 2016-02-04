function cell_type = get_cell_type(datarun, vision_id)

get_cell_indices(datarun, vision_id);

for i = 1:size(datarun.cell_types,2)
    if sum(vision_id == datarun.cell_types{i}.cell_ids)>0
        cell_type = datarun.cell_types{i}.name;
    end
    
end
if ~exist('cell_type')
    cell_type = '';
end

