function [cellIdsToCheck, cellIndices, elecsToCheck] = getLargeAmpSpikes(datarun,spikeThresh)
% Finds cells within a preparation that have spikes with amplitudes below
% specified threshold. returns the cell IDs, indices, and electrodes to
% check. LG

numNeurons         = size(datarun.ei.eis,1); 
allEIs             = zeros(numNeurons, size(datarun.ei.eis{1},1), size(datarun.ei.eis{1},2));
for n  = 1:numNeurons
    allEIs(n,:,:) = datarun.ei.eis{n}; 
end
[minVals,minElecs] = min(min(allEIs,[],3),[],2);
cellIndices        = find(minVals<spikeThresh); 
cellIdsToCheck     = datarun.cell_ids(cellIndices); 
elecsToCheck       = minElecs(cellIndices); 
end