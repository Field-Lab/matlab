function datarun = sync_cell_ids(datarun, cell_ids, source)
% SYNC_CELL_IDS     Ensure that datarun has the correct list of cell IDs
%
% if datarun.cell_ids doesn't match the provided cell_ids, give an error
% if datarun.cell_ids doesn't exist, set it to match the provided cell_ids
%
%
% usage:  datarun = sync_cell_ids(datarun, cell_ids, source)
%
%
% arguments:  datarun - datarun struct with field datarun.cell_ids
%            cell_ids - list of cell IDs
%              source - string specifying the source of cell_ids (e.g. 'STA file')
%
% outputs:    datarun - result of computation
%
%
%
%



% if there is no cell_ids field in datarun, set datarun.cell_ids to match the provided list
if ~isfield(datarun,'cell_ids')
    % ensure it's a row vector
    if size(cell_ids,2) == 1
        cell_ids = cell_ids';
    end
    
    datarun.cell_ids = cell_ids;

else % if there is a datset.cell_ids field, compare it to the provided list

    % if they are not identical
    if length(cell_ids) ~= length(datarun.cell_ids) || ~all(sort(cell_ids) == sort(datarun.cell_ids))

        % give an error
        error('datarun.cell_ids (%d cells) does not match list of cell ids from %s (%d cells).',...
            length(datarun.cell_ids),source,length(cell_ids))
    end
end


