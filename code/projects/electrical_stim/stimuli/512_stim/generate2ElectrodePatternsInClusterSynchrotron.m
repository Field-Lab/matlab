function Array = generate2ElectrodePatternsInClusterSynchrotron(electrodes,sameAmps)

% given a set of a electrodes on the 512 array, this function determines
% which pairs of electrodes within the cluster are nearest-neighbors (30 or
% 60 microns apart for the 30 and 60 micron arrays, respectively), and then
% produces 1 patterns for each of these pairs: 0.5:1
% arguments:
%
%   electrodes, a vector of electrode numbers from the 512 array


nElectrodes = length(electrodes);

if nElectrodes == 7
    nCombinations = 12;
elseif nElectrodes == 5
    nCombinations = 7;
elseif nElectrodes == 4
    nCombinations = 5;
end

combinations = nearestNeighborsFinder512(electrodes);

Array = zeros(nElectrodes, nCombinations*2);
startIndex = 0;
if sameAmps
    for i = 1:nCombinations
        Array(combinations(i,1), 2*(i-1)+1 + startIndex : 2*i + startIndex) = [1 1];
        Array(combinations(i,2), 2*(i-1)+1 + startIndex : 2*i + startIndex) = [1 1];
    end
else
    
    for i = 1:nCombinations
        Array(combinations(i,1), 2*(i-1)+1 + startIndex : 2*i + startIndex) = [0.5 1];
        Array(combinations(i,2), 2*(i-1)+1 + startIndex : 2*i + startIndex) = [1 0.5];
    end
end
