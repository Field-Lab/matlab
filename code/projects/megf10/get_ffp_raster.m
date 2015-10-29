function [raster, raster_all] = get_ffp_raster(datarun, cell_id, duration)

n = length(datarun.triggers)/4;
trigger = zeros(1, n);
for i = 1:n
    trigger(i) = datarun.triggers(4*(i-1)+1);
end
% trigger = trigger - duration/2;

raster = cell(length(cell_id), 1);
raster_all = cell(length(cell_id), 1);
for i = 1:length(cell_id)
    if ~isempty(intersect(datarun.cell_ids, cell_id(i)))
        idx = get_cell_indices(datarun, cell_id(i));
        raster{i, 1} = get_raster(datarun.spikes{idx}, trigger, 'plot', false);
        raster_all{i} = cell2mat(raster{i});
    end
end