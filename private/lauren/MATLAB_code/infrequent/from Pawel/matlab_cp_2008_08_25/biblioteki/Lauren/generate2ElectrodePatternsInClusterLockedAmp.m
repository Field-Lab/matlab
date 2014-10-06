function Array = generate2ElectrodePatternsInClusterLockedAmp(electrodes, individualStim)

% given a cluster of electrodes on the 61 array, this function determines
% which pairs of electrodes within the cluster are nearest-neighbors (30 or
% 60 microns apart for the 30 and 60 micron arrays, respectively), and then
% produces 4 patterns for each of these pairs: 1 1, 1 -1, -1 1, -1 -1.

% arguments 
%   electrodes: a vector of the electrode numbers in the cluster
%   individualStim: binary value specifying whether stimulation with individual electrodes should
%   be included in addition to the paired electrodes (1) or not (0)


nElectrodes = length(electrodes);

% for cases when nLayers = 1
if nElectrodes == 7
    nCombinations = 12;
elseif nElectrodes == 5
    nCombinations = 7;
elseif nElectrodes == 4
    nCombinations = 5;
% for cases when nLayers = 2
elseif nElectrodes == 19
    nCombinations = 42;
elseif nElectrodes == 16
    nCombinations = 34;
elseif nElectrodes == 14
    nCombinations = 29;
elseif nElectrodes == 12
    nCombinations = 23;
elseif nElectrodes == 11
    nCombinations = 21;
elseif nElectrodes == 9
    nCombinations = 16;
else
    error('Unexpected number of electrodes in cluster.  Exiting.')
end

combinations = nearestNeighborsFinder(electrodes, nCombinations);

if individualStim
    Array = zeros(nElectrodes, nCombinations*4 + nElectrodes*2);
else
    Array = zeros(nElectrodes, nCombinations*4);
end


for i = 1:nCombinations
    Array(combinations(i,1), 4*(i-1)+1 : 4*i) = [1 1 -1 -1];
	Array(combinations(i,2), 4*(i-1)+1 : 4*i) = [1 -1 1 -1];
end

if individualStim
    for i = 1:nElectrodes
        Array(i, nCombinations*4 + 2*(i-1) + 1) = 1;
        Array(i, nCombinations*4 + 2*(i-1) + 2) = -1;
    end
end