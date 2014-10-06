function Array = generate2ElectrodePatternsInClusterLockedAmp(electrodes)

% given a set of a electrodes on the 61 array, this function determines
% which pairs of electrodes within the cluster are nearest-neighbors (30 or
% 60 microns apart for the 30 and 60 micron arrays, respectively), and then
% produces 4 patterns for each of these pairs: 1 1, 1 -1, -1 1, -1 -1.

nElectrodes = length(electrodes);

if nElectrodes == 7
    nCombinations = 12;
elseif nElectrodes == 5
    nCombinations = 7;
elseif nElectrodes == 4
    nCombinations = 5;
end

combinations = nearestNeighborsFinder(electrodes, nCombinations);

Array = zeros(nElectrodes, nCombinations*4);
for i = 1:nCombinations
    Array(combinations(i,1), 4*(i-1)+1 : 4*i) = [1 1 -1 -1];
	Array(combinations(i,2), 4*(i-1)+1 : 4*i) = [1 -1 1 -1];
end
