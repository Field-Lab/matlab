function neighbors = hexNeighborsFast(electrodes, hexArray, hexCoords)

neighbors = [];

rmax = size(hexArray, 1);
qmax = size(hexArray, 2);

for electrode = electrodes
    q = hexCoords(electrode, 1);
    r = hexCoords(electrode, 2);
    
    if q-1 >= 1
        neighbor = hexArray(r, q-1);
        if neighbor
            neighbors = cat(2, neighbors, neighbor);
        end
    end
    
    if q+1 <= qmax
        neighbor = hexArray(r, q+1);
        if neighbor
            neighbors = cat(2, neighbors, neighbor);
        end
    end
    
    if r-1 >= 1
        neighbor = hexArray(r-1, q);
        if neighbor
            neighbors = cat(2, neighbors, neighbor);
        end
    end
    
    if r+1 <= rmax
        neighbor = hexArray(r+1, q);
        if neighbor
            neighbors = cat(2, neighbors, neighbor);
        end
    end
    
    if r-1 >= 1 && q+1 <= qmax
        neighbor = hexArray(r-1, q+1);
        if neighbor
            neighbors = cat(2, neighbors, neighbor);
        end
    end
    
    if r+1 <= rmax && q-1 >= 1
        neighbor = hexArray(r+1, q-1);
        if neighbor
            neighbors = cat(2, neighbors, neighbor);
        end
    end
end
neighbors = cat(2, neighbors, electrodes);
neighbors = unique(neighbors);
end