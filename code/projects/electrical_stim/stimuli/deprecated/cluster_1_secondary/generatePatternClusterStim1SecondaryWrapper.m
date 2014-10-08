function [electrodes Array clusterIDs] = generatePatternClusterStim1SecondaryWrapper(centerElectrodes, relAmps, maxDistance)

% generates arrays that will be saved to create 'electrodes' and 'patterns' stimulus files for
% STIM64
%
%
%
% arguments
%    centerElectrodes: vector of electrodes that will be centers of clusters (cannot be on edge of
%    array, and clusters cannot be overlapping)
%    relAmps: a vector of relative secondary electrode amplitudes, not including 0 (primary and
%    secondary alone stimuli added automatically)
%
% returns
%   electrodes, an array of electrode numbers that are used in the
%   generated patterns
%
%   Array, an array of patterns, in which the first dimension corresponds
%   with an electrode in "electrodes" and the second dimension corresponds
%   with the pattern number.




clusterElectrodes = cell(1, length(centerElectrodes));
for i = 1:length(centerElectrodes)
    clusterElectrodes{i} = getCluster(centerElectrodes(i), 'maxDistance', maxDistance);
%     if length(clusterElectrodes{i})<7
%         error('One of the chosen center electrodes is on the edge of the array.  Aborting.')
%     end
end

if(clusterOverlapCheck61(clusterElectrodes))
    error('Clusters overlap.  Aborting.')
end

arrayForEachCluster = cell(1, length(centerElectrodes));

for i = 1:length(centerElectrodes)
    arrayForEachCluster{i} = generatePatternClusterStim1Secondary(centerElectrodes(i),...
        clusterElectrodes{i}, relAmps);
    %arrayForEachCluster{i} = generatePatternClusterStim1Secondary(clusterElectrodes{i});
end

[Array electrodes clusterIDs] = concatenatePatterns(arrayForEachCluster, clusterElectrodes);