centerElectrodes=[32 21];
%no_of_patterns=50;


%coordinate of each electrode on the 61-electrode array, where 1 = 30
%microns for the 60-micron arrays, and 1 = 15 microns for the 30-micron arrays (any edge between
%nearest neighboring electrodes = 2)
%***electrodes 9, 25 and 57 correspond with nans
% originally written by Anastacia, rewritten by Justin Elstrott
% The position of electrode #41 is taken to be the origin (0,0).
% X+ axis: from electrode 41 towards electrode 13.
% Y+ axis: from electrode 41 towards electrode 62.
xCoords = [1.7320508 3.4641016 3.4641016 1.7320508 5.196152 3.4641016 5.196152 6.928203 nan ...
    6.928203 5.196152 3.4641016 6.928203 5.196152 1.7320508 6.928203 3.4641016 5.196152 6.928203 ...
    5.196152 3.4641016 1.7320508 3.4641016 1.7320508 nan 0 1.7320508 0 0 0 -1.7320508 -1.7320508 ...
    -1.7320508 -3.4641016 -3.4641016 -1.7320508 -5.196152 -3.4641016 -5.196152 -6.928203 0 ...
    -6.928203 -5.196152 -3.4641016 -6.928203 -5.196152 -1.7320508 -6.928203 -3.4641016 -5.196152 ...
    -6.928203 -5.196152 -3.4641016 -1.7320508 -3.4641016 -1.7320508 nan 0 -1.7320508 0 0 0 ...
    1.7320508 1.7320508];
yCoords = [3 6 4 1 5 2 3 4 nan 2 1 0 0 -1 -1 -2 -2 -3 -4 -5 -4 -3 -6 -5 nan -2 -7 -4 -6 -8 -7 -5 ...
    -3 -6 -4 -1 -5 -2 -3 -4 0 -2 -1 0 0 1 1 2 2 3 4 5 4 3 6 5 nan 2 7 4 6 8 7 5];


% cd C:\pawel\pliki\nauka\matlab\; 
% fid = fopen('electrodes','wb')
% fwrite(fid,electrodes,'integer*4');
% fclose(fid);

distances = zeros(length(centerElectrodes),64);
for i = 1:length(centerElectrodes)
    for j = 1:64
        distances(i,j) = norm([xCoords(centerElectrodes(i)) - xCoords(j), yCoords(centerElectrodes(i)) - yCoords(j)]);
    end
end

%determines which electrodes lie within the each cluster of 7 (or less if center is on edge of
%array)
clusterElectrodes = cell(1,length(centerElectrodes));
neighborElectrodes = cell(1,length(centerElectrodes));
centerIndex = zeros(1,length(centerElectrodes));
for i = 1:length(centerElectrodes)
    clusterElectrodes{i} = find(squeeze(distances(i,:))<2.1); %includes center electrode and all nearest neighbors
    neighborElectrodes{i} = clusterElectrodes{i};
    centerIndex(i) = find(neighborElectrodes{i}==centerElectrodes(i));
    neighborElectrodes{i}(centerIndex(i)) = [];
end



%checks to make sure clusters don't overlap
for i = 1:64
    for j = 1:length(centerElectrodes)-1
        for k = j+1:length(centerElectrodes)
            if ~isempty(find(i==clusterElectrodes{j},1)) && ~isempty(find(i==clusterElectrodes{k},1))
                error('Clusters overlap.  Aborting.')
            end
        end
    end
end

arrayForEachCluster = cell(length(centerElectrodes));

%generates order of 2-electrode stimulation patterns for each cluster of electrodes
for i = 1:length(centerElectrodes)
    nSurroundingElectrodes = length(clusterElectrodes{i}) - 1;
    arrayForEachCluster{i} = zeros(nSurroundingElectrodes, nSurroundingElectrodes*9);
    arrayForEachCluster{i}(centerIndex(i),:) = ones(1,nSurroundingElectrodes*9);
    electrodeOrderRandomizer = randperm(nSurroundingElectrodes*9);
    electrodeOrder = ceil(electrodeOrderRandomizer/9);
    for j = 1:nSurroundingElectrodes
        relativeAmplitudeOrder = (randperm(9) - 5*ones(1,9))*.25;
        currentElectrode = neighborElectrodes{i}(j);
        currentElectrodeIndexWithinClusterElectrodes = find(clusterElectrodes{i}==currentElectrode);
        electrodeIndecesWithinOrder = find(electrodeOrder==j);
        for k = 1:9
            arrayForEachCluster{i}(currentElectrodeIndexWithinClusterElectrodes, electrodeIndecesWithinOrder(k)) = relativeAmplitudeOrder(k);
        end
    end
end

%interleaves stimulation patterns of the different clusters and shifts electrode indeces
electrodes = [];
maxClusterSize = 0;
for i = 1:length(centerElectrodes)
    electrodes = [electrodes clusterElectrodes{i}]; %#ok<AGROW>
    maxClusterSize = max([maxClusterSize length(clusterElectrodes{i})]);
end

Array = zeros(length(electrodes), (maxClusterSize-1)*length(centerElectrodes)*9);

placeHolder = 0;
for i = 1:length(centerElectrodes)
    for j = 1:size(arrayForEachCluster{i},2)
        currentArrayPatternIndex = length(centerElectrodes)*(j-1) + i;
        Array(placeHolder+1 : placeHolder+size(arrayForEachCluster{i},1), currentArrayPatternIndex) = arrayForEachCluster{i}(:,j);
    end
    placeHolder = placeHolder + size(arrayForEachCluster{i},1);
end


%Array=rand(length(electrodes),no_of_patterns);
cd C:\pawel\pliki\nauka\matlab\; 
fid = fopen('patterns_Lauren','wb','ieee-le.l64')
fwrite(fid,Array,'double');
fclose(fid);