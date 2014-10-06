function clusterElectrodes = getCluster(centerElectrode, nLayers)

% determines which electrodes lie within each cluster of 7 (nLayers = 1) or 42 (nLayers = 2) (or 
% less if center is on/near edge of array)
%
% puts centerElectrode first in the order
%
% arguments
%   centerElectrode: scalar denoting which electrode on the 61-electrode array is the center of the
%   cluster
%   nLayers: number of "layers" of electrodes surrounding each center that should be included in the
%   clusters (1 gives 4-7 electrodes in cluster, 2 gives 19 electrodes in cluster) (default = 1)

if nargin < 2
    nLayers = 1;
end

[xCoords yCoords] = getElectrodeCoords61();

distances = zeros(1, 64);

for i = 1:64
    distances(i) = norm([xCoords(centerElectrode) - xCoords(i), yCoords(centerElectrode) - yCoords(i)]);
end

if nLayers == 1
    clusterElectrodesTemp = find(squeeze(distances)<2.1); %includes center electrode and all nearest neighbors
else
    clusterElectrodesTemp = find(squeeze(distances)<4.1);
end
neighborElectrodes = clusterElectrodesTemp;
neighborElectrodes(neighborElectrodes==centerElectrode) = [];
clusterElectrodes = [centerElectrode neighborElectrodes]; %puts center electrode first in the order

