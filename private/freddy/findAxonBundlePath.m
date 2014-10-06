function path = findAxonBundlePath(voltageMap, borderElecs)

voltageMap = abs(voltageMap);

gscore = zeros(512, 1);
fscore = zeros(512, 1);

%find the electrode with the greatest voltage (abs) along the edge
maxElec = 0;
maxVoltage = 0;

hexCoords = elecHexCoords();

topBorder = [249:256 261:8:381 385:392];
rightBorder = 392:8:512;
bottomBorder = [505:512 8:8:136 129:135];
leftBorder = 129:8:249;
borderElecs = unique([topBorder rightBorder bottomBorder leftBorder]);
borderVoltages = voltageMap(borderElecs);

[~, ind] = max(borderVoltages);
maxElec = borderElecs(ind);

excludeNeighbors = maxElec;
for i = [1, 2]
    excludeNeighbors = union(excludeNeighbors, intersect(hexNeighborsForElec(excludeNeighbors, hexCoords), borderElecs));
end
borderElecs(ismember(borderElecs, excludeNeighbors)) = [];
borderVoltages = voltageMap(borderElecs);

[~, ind] = max(borderVoltages);
destElec = borderElecs(ind);

% if ismember(maxElec, topBorder)
%     destElec = find(voltageMap == max(voltageMap(bottomBorder)), 1);
% elseif ismember(maxElec, rightBorder)
%     destElec = find(voltageMap == max(voltageMap(leftBorder)), 1);
% elseif ismember(maxElec, bottomBorder)
%     destElec = find(voltageMap == max(voltageMap(topBorder)), 1);
% elseif ismember(maxElec, leftBorder)
%     destElec = find(voltageMap == max(voltageMap(rightBorder)), 1);
% end

maxVoltage = max(voltageMap, [], 2);

weightedMap = (-voltageMap + maxVoltage)/maxVoltage;

destQ = hexCoords(destElec, 1);
destR = hexCoords(destElec, 2);

visited = zeros(512, 1);
openset = zeros(512, 1);
camefrom = zeros(512,1);

gscore(maxElec) = 0;
fscore(maxElec) = gscore(maxElec) + elecHexDist(hexCoords(maxElec, 1),hexCoords(maxElec, 2), destQ, destR)*1.001;

openset(maxElec) = 1;

while any(openset)
    current = find(fscore == min(fscore(openset == 1)), 1);
    if current == destElec
        path = reconstructPath(camefrom, destElec);
        break;
    end
    
    openset(current) = 0;
    visited(current) = 1;
    neighbors = hexNeighborsForElec(current, hexCoords);
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
    
end

    
    
    
    