% Analyze pulse triplet experiments in rat data 2014-07-24-0
piece = '/Volumes/Analysis/2014-07-24-0/';

datarun = load_data([piece 'data000-lg/data000-lg']);
datarun = load_params(datarun,'verbose',1);

% datarun = load_sta(datarun,'verbose',1,'load_sta',{1});
% datarun = get_sta_summaries(datarun,{1});
% cellSpec = get_cell_indices( datarun, {'ON parasol'} );
cell_specification = datarun.cell_ids(:); %Trying first cell
datarun = load_ei(datarun, 'all');


% Fine cells with large spikes
spikeThresh = -80; 
eiFile = edu.ucsc.neurobiology.vision.io.PhysiologicalImagingFile(datarun.names.rrs_ei_path);
numNeurons = size(datarun.ei.eis,1); 
allEIs = zeros(numNeurons, size(datarun.ei.eis{1},1), size(datarun.ei.eis{1},2));
for i  = 1:numNeurons
    allEIs(i,:,:) = datarun.ei.eis{i}; 
end

mins        = min(min(allEIs,[],3),[],2);
cellIndices = find(mins<spikeThresh); 
cellsToCheck = cell_specification(cellIndices); 

temp = min(allEIs(cellIndices,:,:),[],3);
[~,elecsToCheck] = find(temp == repmat(min(temp,[],2),1,512));


