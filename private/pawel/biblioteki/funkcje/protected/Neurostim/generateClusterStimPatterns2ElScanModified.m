function [electrodes Array] = generateClusterStimPatterns2ElScanModified(centerElectrodes)

clusterElectrodes = cell(1, length(centerElectrodes));
for i = 1:length(centerElectrodes)
    clusterElectrodes{i} = getCluster(centerElectrodes(i));
    if length(clusterElectrodes{i})<7
        error('One of the chosen center electrodes is on the edge of the array.  Aborting.')
    end
end

if(clusterOverlapCheck61(clusterElectrodes))
    error('Clusters overlap.  Aborting.')
end

arrayForEachCluster = cell(1, length(centerElectrodes));

for i = 1:length(centerElectrodes)
    arrayForEachCluster{i} = generate2ElectrodePatternsInClusterModified(clusterElectrodes{i});
end

[Array electrodes] = concatenatePatterns(arrayForEachCluster, clusterElectrodes);