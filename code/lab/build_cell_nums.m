function datarun = build_cell_nums(datarun)
% BUILD_CELL_NUMS    Create an index table for reverse lookup of cell_num by cell_id
%
% Simply inverts the lookup throught the cell_id table.  Uses sparse matrix
% for efficiency, although this is probably not necessary.
%
% 2010-01 phli

datarun.cell_nums = sparse([], [], [], double(max(datarun.cell_ids)), 1, numel(datarun.cell_ids));

for i = 1:numel(datarun.cell_ids)
    datarun.cell_nums(datarun.cell_ids(i)) = i;
end