function summ = the_rf(datarun, cell_id)
% SUMMARY     get summary frame from specified cell ID
%
% usage:  datarun = summary(datarun, cell_id)
%
% arguments:  datarun - datarun struct with field datarun.stas.rfs
%             cell_id - cell id
%
% outputs:       summ - summary frame of the desired cell
%
%


% BODY OF THE FUNCTION

summ = datarun.stas.rfs{get_cell_indices(datarun,cell_id)};

