function [X, Y] = getCellEllipse(datarun, neuronID)

[cell_index,~,~] = get_cell_indices(datarun, neuronID);

the_fit = datarun.stas.fits{cell_index};
ctr = the_fit.mean;
rad = the_fit.sd;
[X,Y] = drawEllipse([ctr rad the_fit.angle]);

end