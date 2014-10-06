function clusterElectrodes = getCluster(centerElectrode)

% determines which electrodes lie within the each cluster of 7 (or less if center is on edge of
% array)
% puts centerElectrode first in the order

[xCoords yCoords] = getElectrodeCoords61();

distances = zeros(1, 64);

for i = 1:64
    distances(i) = norm([xCoords(centerElectrode) - xCoords(i), yCoords(centerElectrode) - yCoords(i)]);
end

clusterElectrodesTemp = find(squeeze(distances)<2.1); %includes center electrode and all nearest neighbors
neighborElectrodes = clusterElectrodesTemp;
neighborElectrodes(find(neighborElectrodes==centerElectrode)) = [];
clusterElectrodes = [centerElectrode neighborElectrodes]; %puts center electrode first in the order

