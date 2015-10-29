function raster_all = get_mb_raster(datarun, cell_id, duration)
% raster_all = get_ds_raster(datarun, cell_id)
% datarun: drifting grating dataset
% raster_all: S cells, M spatial period x N temporal period x R direction x P repeats
% for each cell
%
% xyao
% 2013-12-16

triggers = datarun.stimulus.triggers;
index = grp_stim_mb(datarun);

raster_all = cell(length(cell_id), 1);
for rgc = 1:length(cell_id)
    if ismember(cell_id(rgc), datarun.cell_ids)
    idx = get_cell_indices(datarun, cell_id(rgc));
    spike = datarun.spikes{idx};
    raster = get_raster_xy(spike, triggers, duration, 'plot', false);
    raster_all{rgc} = raster(index);
    end
end

    
                

