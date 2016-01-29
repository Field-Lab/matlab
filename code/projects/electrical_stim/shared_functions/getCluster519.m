function clusterElectrodes = getCluster519(centerElectrode, varargin)

% determines which electrodes lie within the each cluster of 7 (or less if center is on edge of
% array)
% puts centerElectrode first in the order
% LG modified code from LHJ to analyze data from the 512 electrode array.
% 1/24/2014

p = inputParser;
p.addRequired('centerElectrode', @isnumeric)
p.addParamValue('maxDistance', 1, @isnumeric) %default: nearest-neighbors only
p.parse(centerElectrode, varargin{:})

maxDist = p.Results.maxDistance*30; %factor of 60 because closest neighbor distance in getElectrodeCoords512 is 2

%%
[xCoords, yCoords] = getElectrodeCoords519();
distances = zeros(1, size(xCoords,2));

for i = 1:size(xCoords,2)
    distances(i) = norm([xCoords(centerElectrode) - xCoords(i), yCoords(centerElectrode) - yCoords(i)]);
end

clusterElectrodesTemp = find(squeeze(distances) < maxDist*1.15); %includes center electrode and all nearest neighbors
neighborElectrodes = clusterElectrodesTemp;
neighborElectrodes(neighborElectrodes==centerElectrode) = [];
clusterElectrodes = [centerElectrode neighborElectrodes]; %puts center electrode first in the order

