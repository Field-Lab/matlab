function datarun = make_unclassified_cell_type(datarun)
% make_unclassified_cell_type     make a new cell class for cells that are not classified
%
% usage:  datarun = make_unclassified_cell_type(datarun)
%
% arguments:  datarun - datarun struct with field datarun.cell_ids
%
% outputs:    datarun - datarun struct
%
% 2008-11 gauthier
%


% make an extra class for cells that are not classified

%accumulate list of classified cells
classified_cells = [];
for cc =1:length(datarun.cell_types)
    classified_cells = [classified_cells datarun.cell_types{cc}.cell_ids];

end

% identify unclassified cells
leftover_cells = setdiff(datarun.cell_ids,classified_cells);

% if any put them in a new class
if ~isempty(leftover_cells)
    class_index = length(datarun.cell_types)+1;
    datarun.cell_types{class_index}.cell_ids = leftover_cells;
    datarun.cell_types{class_index}.name = 'unclassified';
end

