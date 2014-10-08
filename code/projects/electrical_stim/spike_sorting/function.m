function neuronsAboveThresh = findNeuronsAboveThresh(eiPath, channels, thresh, neurons)

if ~exist('neurons', 'var')
    neurons = 1;
end

nNeurons = length(neurons);
neuronMaxs = zeros(nNeurons, 1);

eiFile = edu.ucsc.neurobiology.vision.io.PhysiologicalImagingFile(eiPath);
for i = 1:nNeurons
    eiData = eiFile.getImage(neurons(i));
    eiData = squeeze(eiData{i}(1, 2:end, :));
    
    neuronMaxs(i) = max(max(abs(eiData(channels, :))));
end

neuronsAboveThresh = neurons(neuronMaxs >= thresh);