function [electrodes Array] = generatePatternsAxonStim(centerElectrodes, individualStim)

% arguments
%   centerElectrodes: a vector of electrodes on the 61 electrode array
%   individualStim: binary value specifying whether stimulation with individual electrodes should
%   be included in addition to the paired electrodes (1) or not (0)
%
%
% returns: 
%
%   electrodes, an array of electrode numbers that are used in the
%   generated patterns
%
%   Array, an array of patterns, in which the first dimension corresponds
%   with an electrode in "electrodes" and the second dimension corresponds
%   with the pattern number.
%
%   see help for generate2ElectrodePatternsInClusterLockedAmp for details
%   about patterns generated

nClusters = length(centerElectrodes);

% if nLayers == 1
    electrodes = zeros(7*nClusters, 1);
    Array = zeros(7*nClusters, 62*nClusters);
    clusterElectrodes = cell(nClusters);
    clusterArrays = cell(nClusters);
    for i = 1:nClusters
        %generates electrodes and Arrays
        clusterElectrodes{i} = getCluster(centerElectrodes(i));
        clusterArrays{i} = generate2ElectrodePatternsInClusterLockedAmp(clusterElectrodes{i}, individualStim);
        %concatenates electrodes and Arrays; maximum # of patterns for one cluster = 62 (includes
        %individuals, cluster of 7 electrodes)
        electrodes(7*(i-1)+1 : 7*(i-1) + length(clusterElectrodes{i})) = clusterElectrodes{i};
        Array(7*(i-1)+1 : 7*(i-1) + size(clusterArrays{i}, 1), 62*(i-1) + 1 : 62*(i-1) + size(clusterArrays{i}, 2)) = clusterArrays{i};
    end
    if(clusterOverlapCheck61(clusterElectrodes))
        error('Clusters overlap.  Aborting.')
    end
% else
%     electrodes = getCluster(centerElectrodes, nLayers);
%     Array = generate2ElectrodePatternsInClusterLockedAmp(electrodes, individualStim);
% end