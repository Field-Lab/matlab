

centerElectrodes=[32 21];


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


cd C:\pawel\pliki\nauka\matlab\; 
fid = fopen('electrodes','wb') %#ok<NOPTS>
fwrite(fid,electrodes,'integer*4');
fclose(fid);


distances = zeros(length(centerElectrodes), 64);
for i = 1:length(centerElectrodes)
    for j = 1:64
        distances(i,j) = norm([xCoords(centerElectrodes(i)) - xCoords(j), yCoords(centerElectrodes(i)) - yCoords(j)]);
    end
end

%determines which electrodes lie within the each cluster of 7 (or less if center is on edge of
%array)
clusterElectrodes = cell(1,length(centerElectrodes));
neighborElectrodes = cell(1,length(centerElectrodes));
for i = 1:length(centerElectrodes)
    clusterElectrodesTemp = find(squeeze(distances(i,:))<2.1); %includes center electrode and all nearest neighbors
    neighborElectrodes{i} = clusterElectrodesTemp;
    centerIndex = find(neighborElectrodes{i}==centerElectrodes(i));
    neighborElectrodes{i}(centerIndex) = [];
    clusterElectrodes{i} = [centerElectrodes(i) neighborElectrodes{i}]; %puts center electrode first in the order
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

for i = 1:length(centerElectrodes)
    nSurroundingElectrodes = length(neighborElectrodes{i});
    arrayForEachCluster{i} = zeros(length(clusterElectrodes{i}), nSurroundingElectrodes*16+2);
    arrayForEachCluster{i}(1,1) = 1;
    arrayForEachCluster{i}(1,2) = -1;
    for j = 1:4
      for k = 1:nSurroundingElectrodes
          arrayForEachCluster{i}(1+k, 2+k+(j-1)*nSurroundingElectrodes*4) = 0.25*j;
          arrayForEachCluster{i}(1, 2+k+(j-1)*nSurroundingElectrodes*4) = 1;
          arrayForEachCluster{i}(1+k, 2+k+(j-1)*nSurroundingElectrodes*4 + nSurroundingElectrodes) = 0.25*j;
          arrayForEachCluster{i}(1, 2+k+(j-1)*nSurroundingElectrodes*4 + nSurroundingElectrodes) = -1;
          arrayForEachCluster{i}(1+k, 2+k+(j-1)*nSurroundingElectrodes*4 + nSurroundingElectrodes*2) = -0.25*j;
          arrayForEachCluster{i}(1, 2+k+(j-1)*nSurroundingElectrodes*4 + nSurroundingElectrodes*2) = 1;
          arrayForEachCluster{i}(1+k, 2+k+(j-1)*nSurroundingElectrodes*4 + nSurroundingElectrodes*3) = -0.25*j;
          arrayForEachCluster{i}(1, 2+k+(j-1)*nSurroundingElectrodes*4 + nSurroundingElectrodes*3) = -1;
      end
    end
end


electrodes = [];
Array = [];
for i = 1:length(centerElectrodes)
    electrodes = [electrodes clusterElectrodes{i}]; %#ok<AGROW>
    Array = [Array arrayForEachCluster{i}]; %#ok<AGROW>
end


fid = fopen('patterns','wb','ieee-le.l64') %#ok<NOPTS>
fwrite(fid,Array,'double');
fclose(fid);