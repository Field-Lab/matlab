function [electrodes Array] = generateClusterStimPatterns2ElScan(centerElectrodes)


clusterElectrodes = cell(1, length(centerElectrodes));
for i = 1:length(centerElectrodes)
    clusterElectrodes{i} = getCluster(centerElectrodes(i));
end

if(clusterOverlapCheck61(clusterElectrodes))
    error('Clusters overlap.  Aborting.')
end

arrayForEachCluster = cell(1, length(centerElectrodes));

for i = 1:length(centerElectrodes)
    arrayForEachCluster{i} = generate2ElectrodePatternsInCluster(clusterElectrodes{i});
end

[Array electrodes] = concatenatePatterns(arrayForEachCluster, clusterElectrodes);