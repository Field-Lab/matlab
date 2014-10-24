function init = mergeRows(init,rows)

% merge two rows of an assignment matrix

tmp = init(rows(1),:) | init(rows(2),:);
init(rows(1),:) = tmp;
init(rows(2),:) = tmp;
init = removeSubsets(init);