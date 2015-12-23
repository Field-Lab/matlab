function cell_nums = nums_from_ids(datarun, cell_ids)
% convert cell ids to row vector of cell numbers
%
% 2010-01 phli
%


cell_nums = zeros(1,length(cell_ids));
for cc = 1:length(cell_ids)

    % look for the specified cell id in the list of cell ids
    temp = find(datarun.cell_ids == cell_ids(cc));

    if ~isempty(temp)
        % if it was found, add it to the list of cell numbers
        cell_nums(cc) = temp;
    else
        % if it wasn't found, give an error
        error('cell id %d does not exist',cell_ids(cc))
    end
end