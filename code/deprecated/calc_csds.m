function datarun = calc_csds(datarun, cellspec)
if ~isfield(datarun.ei, 'neighbor_struct')
    datarun = get_lazyhex_ei_neighbors(datarun);
end
neighbor_struct = datarun.ei.neighbor_struct;
cellnums = get_cell_indices(datarun, cellspec);
datarun.ei.csd(cellnums) = cellfun(@(ei)(ei2csd(ei,neighbor_struct)), datarun.ei.eis(cellnums), 'UniformOutput', false);