function pols = get_polarities(datarun, cellspec)
cellnums = get_cell_indices(datarun, cellspec);
pols = [datarun.stas.polarities{cellnums}];