function [Array electrodes] = concatenatePatterns(arrayForEachCluster, clusterElectrodes)

electrodes = [clusterElectrodes{1}];
Array = arrayForEachCluster{1};

if length(clusterElectrodes) > 1
    for i = 2:length(clusterElectrodes)
        electrodes = [electrodes clusterElectrodes{i}]; %#ok<AGROW>
        spacer1 = zeros(size(arrayForEachCluster{i}, 1), size(Array, 2));
        spacer2 = zeros(size(Array, 1), size(arrayForEachCluster{i}, 2));
        Array = cat(1, Array, spacer1);
        toBeCatted = cat(1, spacer2, arrayForEachCluster{i});
        Array = cat(2, Array, toBeCatted);
    end
end