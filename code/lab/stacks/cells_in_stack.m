function [cell_ids, cell_nums] = cells_in_stack(datarun, stack_index, cell_method)
% CELLS_IN_STACK    Determine which cells are within stack image bounds
% usage: [cell_ids, cell_nums] = cells_in_stack(datarun, stack_index, cell_method)
% 
% For now, the only cell_method is 'channels', meaning the coordinates of
% each cell are based only on the electrode the cell was detected on.
%
% 2010-08 phli
%

if nargin < 3
    cell_method = 'channels'; % For now; once better options are coded probably should change this default.
end

switch cell_method
    case 'channels'
        primary_electrodes_per_cell = datarun.channels;
        positions = datarun.ei.position(primary_electrodes_per_cell, 1:2);
end

stack = get_stack(datarun, stack_index);
cell_nums = points_in_stack(stack, positions);
cell_ids = datarun.cell_ids(cell_nums);