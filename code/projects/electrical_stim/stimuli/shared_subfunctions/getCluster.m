function clusterElectrodes = getCluster(centerElectrode, varargin)

% determines which electrodes lie within the each cluster of 7 (or less if center is on edge of
% array)
% puts centerElectrode first in the order

p = inputParser;

p.addRequired('centerElectrode', @isnumeric)

p.addParamValue('maxDistance', 1, @isnumeric) %default: nearest-neighbors only

p.parse(centerElectrode, varargin{:})

maxDist = p.Results.maxDistance*2; %factor of 2 because closest neighbor distance in getElectrodeCoords is 2

%%
[xCoords yCoords] = getElectrodeCoords61();

distances = zeros(1, 64);

for i = 1:64
    distances(i) = norm([xCoords(centerElectrode) - xCoords(i), yCoords(centerElectrode) - yCoords(i)]);
end

clusterElectrodesTemp = find(squeeze(distances) < maxDist*1.05); %includes center electrode and all nearest neighbors
neighborElectrodes = clusterElectrodesTemp;
neighborElectrodes(neighborElectrodes==centerElectrode) = [];
clusterElectrodes = [centerElectrode neighborElectrodes]; %puts center electrode first in the order

