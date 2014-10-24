function inds = getCellInds(rgc_ids,ids)

% input: vector of all cell ids, set of ids to get indices for
% output: vector of indices for the desired cells

inds = zeros(size(ids));

for iCell=1:length(ids)
    inds(iCell) = find(ids(iCell)==rgc_ids);
end


