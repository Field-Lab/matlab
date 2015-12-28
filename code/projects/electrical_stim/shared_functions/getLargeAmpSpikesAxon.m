function [cellIdsToCheck, cellIndices, elecsToCheck] = getLargeAmpSpikesAxon(datarun,spikeThresh,varargin)
% Finds cells within a preparation that have spikes with amplitudes below
% specified threshold. returns the cell IDs, indices, and electrodes to
% check. LG

%reformat EI data
numNeurons         = size(datarun.ei.eis,1); 
allEIs             = zeros(numNeurons, size(datarun.ei.eis{1},1), size(datarun.ei.eis{1},2));
for n  = 1:numNeurons
    allEIs(n,:,:) = datarun.ei.eis{n}; 
end

%determine which electrodes somas pass threshold (to avoid junk)
if length(varargin) == 0
    soma_thr = -20; %<---PARAMETER
else 
    soma_thr = varargin{1}
end
[minVals,minElecs] = min(min(allEIs,[],3),[],2);
cellIndices        = find(minVals<soma_thr); 
cellIdsToCheck     = datarun.cell_ids(cellIndices); 
elecsToCheck       = minElecs(cellIndices); 

%of cells which pass threshold, examine axons
load('adj_mat_512.mat');
minMat = min(allEIs,[],3); 
minMat = minMat(cellIndices, :);
for n = 1:length(cellIndices)
    mainElec = elecsToCheck(n); %assuming minimum electrode for this neuron is at soma
    %zero out main electrode, all adjacents, and all adjacents to each adjacent
    minMat(n,mainElec) = 0; 
    for i = adj_mat_512{mainElec}
	minMat(n,i) = 0;
	for j = adj_mat_512{i}
	    minMat(n,j) = 0;
	end
    end
end
[minVals,minElecs] = min(minMat,[],2);
indices        = find(minVals<spikeThresh); %this is now the axonal threshold
cellIdsToCheck     = cellIdsToCheck(indices); 
elecsToCheck       = minElecs(indices); 
cellIndices = indices;
end
