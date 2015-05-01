function path = findAxonBundlePathFast(voltageMap, borderElecs, hexCoords, startElec, endElec, exclude)

voltageMap = abs(voltageMap);

gscore = zeros(512, 1);
fscore = zeros(512, 1);

hexArray = elecHexArray(hexCoords);

% borderElecs(borderElecs == exclude) = [];

%find the electrode with the greatest voltage (abs) along the edge
if nargin < 5
    if nargin < 4
        startElec = 0;
    end
    endElec = 0;
end

if ~startElec
    borderVoltages = voltageMap(borderElecs);

    [~, ind] = max(borderVoltages); 
    maxElec = borderElecs(ind);
    startElec = maxElec;
end

if ~endElec
    if ismember(startElec, borderElecs)
        excludeNeighbors = startElec;
        for i = 1:3
            excludeNeighbors = union(excludeNeighbors, intersect(hexNeighborsFast(excludeNeighbors, hexArray, hexCoords), borderElecs));
        end
        borderElecs(ismember(borderElecs, excludeNeighbors)) = [];
    end
        borderVoltages = voltageMap(borderElecs);
        [~, ind] = max(borderVoltages);
        destElec = borderElecs(ind);
        endElec = destElec;
end

maxVoltage = max(voltageMap, [], 2);

weightedMap = (-voltageMap + maxVoltage)/maxVoltage;

destQ = hexCoords(endElec, 1);
destR = hexCoords(endElec, 2);

visited = zeros(512, 1);
openset = zeros(512, 1);
camefrom = zeros(512,1);

gscore(startElec) = 0;
fscore(startElec) = gscore(startElec) + elecHexDist(hexCoords(startElec, 1),hexCoords(startElec, 2), destQ, destR)*1.001;

openset(startElec) = 1;

while any(openset)
    current = find(fscore == min(fscore(openset == 1)), 1);
    if current == endElec
        path = reconstructPath(camefrom, endElec);
        break;
    end
    
    openset(current) = 0;
    visited(current) = 1;
    neighbors = hexNeighborsFast(current, hexArray, hexCoords);
    for neighbor = neighbors
        if (visited(neighbor))
            continue;
        end
        gestimate = gscore(current) + weightedMap(neighbor);%elecHexDist(hexCoords(current, 1),hexCoords(current, 2), hexCoords(neighbor, 1), hexCoords(neighbor, 2));
        if ~openset(neighbor) || gestimate < gscore(neighbor)
            camefrom(neighbor) = current;
            gscore(neighbor) = gestimate;
            hueristic = elecHexDist(hexCoords(neighbor, 1),hexCoords(neighbor, 2), destQ, destR);
            fscore(neighbor) = gscore(neighbor) + hueristic*1.001;
            openset(neighbor) = 1;
        end
    end
end