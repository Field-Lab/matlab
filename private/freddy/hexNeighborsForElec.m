function neighbors = hexNeighborsForElec(electrodes, hexCoords)
neighbors = [];
for electrode = electrodes
    q = hexCoords(electrode, 1);
    r = hexCoords(electrode, 2);
    neighbor = find(ismember(hexCoords, [q+1, r], 'rows'));
    if neighbor
        neighbors = cat(2, neighbors, neighbor);
    end
    neighbor = find(ismember(hexCoords, [q-1, r], 'rows'));
    if neighbor
        neighbors = cat(2, neighbors, neighbor);
    end
    neighbor = find(ismember(hexCoords, [q, r+1], 'rows'));
    if neighbor
        neighbors = cat(2, neighbors, neighbor);
    end
    neighbor = find(ismember(hexCoords, [q, r-1], 'rows'));
    if neighbor
        neighbors = cat(2, neighbors, neighbor);
    end
    neighbor = find(ismember(hexCoords, [q-1, r+1], 'rows'));
    if neighbor
        neighbors = cat(2, neighbors, neighbor);
    end
    neighbor = find(ismember(hexCoords, [q+1, r-1], 'rows'));
    if neighbor
        neighbors = cat(2, neighbors, neighbor);
    end
end
neighbors = unique(neighbors);
end