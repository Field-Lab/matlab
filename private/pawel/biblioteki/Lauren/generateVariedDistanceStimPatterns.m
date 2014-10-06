function [electrodes Array] = generateVariedDistanceStimPatterns(targetElectrodes, individualStim)

% arguments
%   targetElectrodes: a vector of electrodes on the 61 electrode array
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

nTargets = length(targetElectrodes);

electrodes = (1:64)';

Array = zeros(64, 61+60*nTargets);
%Array = zeros(7*nClusters, 62*nClusters);

%single-electrode stimulation
patternCount = 0;
for i = 1:64
    if (i~=9)&&(i~=25)&&(i~=57)
        patternCount = patternCount+1;
        Array(i,patternCount) = 1;
    end
end

if patternCount~=61
    error('something went wrong with creating individual electrode stim patterns')
end

%paired-electrode stimulation
for i = 1:nTargets
    for j = 1:64
        if (j~=9)&&(j~=25)&&(j~=57)&&(j~=targetElectrodes(i))
            patternCount = patternCount + 1;
            Array(j,patternCount) = 1;
            Array(targetElectrodes(i),patternCount) = 1;
        end
    end
end

if patternCount ~= 61+60*nTargets
    error('something went wrong with creating paired electrode stim patterns')
end






%clusterElectrodes = cell(nClusters);
%clusterArrays = cell(nClusters);




%for i = 1:nClusters
    %generates electrodes and Arrays
    %clusterElectrodes{i} = getCluster(centerElectrodes(i), nLayers);
    %clusterArrays{i} = generate2ElectrodePatternsInClusterLockedAmp(clusterElectrodes{i}, individualStim);
    %concatenates electrodes and Arrays
    %electrodes(7*(i-1) + 1 : 7*(i-1) + length(clusterElectrodes{i})) = clusterElectrodes{i};
    %Array(7*(i-1) + 1 : 7*(i-1) + size(clusterArrays{i}, 1), 62*(i-1) + 1 : 62*(i-1) + size(clusterArrays{i}, 2)) = clusterArrays{i};
%end