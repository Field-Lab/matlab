function show_cell_types(cell_types)
% show_cell_types     Display a list of how many cells are in each type
%
% usage:  show_cell_types(cell_types)
%         show_cell_types(datarun.cell_types)
%
% arguments:  cell_types - standard cell_types struct
%



for tt = 1:length(cell_types)
    %fprintf('\n\t%d\t%s',length(cell_types{tt}.cell_ids),cell_types{tt}.name)
    %fprintf('\n\t%d.  %s (%d cells)',tt,cell_types{tt}.name,length(cell_types{tt}.cell_ids))
    fprintf('\n\t%d\t%s (%d)',length(cell_types{tt}.cell_ids),cell_types{tt}.name,tt)
end
fprintf('\n')

