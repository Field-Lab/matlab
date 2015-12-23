function overlapping = findOverlappingNeurons(datarun, threshold, neurons)

topBorder = [249:256 261:8:381 385:392];
rightBorder = 392:8:512;
bottomBorder = [505:512 8:8:136 129:135];
leftBorder = 129:8:249;
borderElecs = unique([topBorder rightBorder bottomBorder leftBorder]);

hexCoords = elecHexCoords;

if nargin < 3
    [neurons, ~, ~] = getLargeAmpSpikes(datarun, threshold);
end

overlappingLarge = zeros(size(neurons,2), size(neurons,2));

paths = zeros(45, size(neurons,2));
i = 1;
for neuron = neurons
    EI = datarun.ei.eis{get_cell_indices(datarun, neuron)};
    eiamps = max(EI') - min(EI');
    [~, maxIndex] = max(eiamps);
    path = findAxonBundlePathFast(eiamps, borderElecs, hexCoords, maxIndex);
    paths(1:size(path, 2), i) = path;
    i = i + 1;
end

maxOverlaps = 1;
for i = 1:size(neurons,2)
    path1 = paths(:, i);
    overlappingLarge(i, 1) = neurons(i);
    path1(path1 == 0) = [];
    numOverlaps = 1;
    disp(['Checking for overlaps on neuron: ' num2str(neurons(i))]);
    for j= 1:size(neurons,2)
        if neurons(i) == neurons(j)
            continue;
        end
        path2 = paths(:,j);
        path2(path2 == 0) = [];
        if size(path1, 1) > size(path2, 1)
            overlaps = bundlesOverlap(path2, path1);
        else
            overlaps = bundlesOverlap(path1, path2);
        end
        if overlaps
            numOverlaps = numOverlaps + 1;
            overlappingLarge(i, find(overlappingLarge(i, :) == 0, 1, 'first')) = neurons(j);
        end
    end
    overlappingLarge(i, size(overlappingLarge, 2)) = numOverlaps;
    if numOverlaps > maxOverlaps
        maxOverlaps = numOverlaps;
    end
end
overlapping = overlappingLarge(:, 1:(maxOverlaps+1));
overlapping(:, maxOverlaps+1) = overlappingLarge(:, size(overlappingLarge,2));
            
end