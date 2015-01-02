function [electrodes Array] = generatePatternsSynchrotronStim(centerElectrodes, sameAmps)

% arguments
%   centerElectrodes: a vector of electrodes on the 512 electrode array
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

nClusters  = length(centerElectrodes);
electrodes = zeros(7*nClusters, 1);
% Array size is numElectrodes, numPatterns
numAmpRatios = 1; 
numPairs = 12; 
numDirections = 2;
numAmplitudesTotal =2;  % Planning to test a current ratio of 0.5:1
maxNumPatterns = numAmpRatios*numPairs*numDirections + numAmplitudesTotal*7; 

Array = zeros(size(electrodes,1), maxNumPatterns*nClusters);
clusterElectrodes = cell(nClusters);
clusterArrays = cell(nClusters);
for i = 1:nClusters
    %generates electrodes and Arrays
    clusterElectrodes{i} = getCluster512(centerElectrodes(i));
    
    clusterArrays{i} = generate2ElectrodePatternsInClusterSynchrotron(clusterElectrodes{i},sameAmps);
    %concatenates electrodes and Arrays; maximum # of patterns for one cluster = maxNumPatterns (includes
    %cluster of 7 electrodes, no individuals)
    electrodes(7*(i-1)+1 : 7*(i-1) + length(clusterElectrodes{i})) = clusterElectrodes{i};
    Array(7*(i-1)+1 : 7*(i-1) + size(clusterArrays{i}, 1), maxNumPatterns*(i-1) + 1 : maxNumPatterns*(i-1) + size(clusterArrays{i}, 2)) = clusterArrays{i};
    if sameAmps
        Array(7*(i-1)+1 : 7*(i-1) + size(clusterArrays{i}, 1), maxNumPatterns*(i-1) + size(clusterArrays{i}, 2) + 1 : 1 : maxNumPatterns*(i-1) + size(clusterArrays{i}, 2) + 7) = 1*diag(ones(7,1));
        
    else
        Array(7*(i-1)+1 : 7*(i-1) + size(clusterArrays{i}, 1), maxNumPatterns*(i-1) + size(clusterArrays{i}, 2) + 1 : 1 : maxNumPatterns*(i-1) + size(clusterArrays{i}, 2) + 7) = 0.5*diag(ones(7,1));
    end
    Array(7*(i-1)+1 : 7*(i-1) + size(clusterArrays{i}, 1), maxNumPatterns*(i-1) + size(clusterArrays{i}, 2) + 8 : 1 : maxNumPatterns*(i-1) + size(clusterArrays{i}, 2) + 14) = 1*diag(ones(7,1));
end
if(clusterOverlapCheck512(clusterElectrodes))
    error('Clusters overlap.  Aborting.')
end
