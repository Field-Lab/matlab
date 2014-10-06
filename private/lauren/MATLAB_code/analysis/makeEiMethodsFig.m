pathToEi = '/snle/lab/Experiments/Array/Analysis/2008-08-27-2/data001-lh/data001-lh.ei';
neuronID = 811;

figure
axesH = gca;
plotEi61(pathToEi, neuronID, 'axesH', axesH)




eiFile = edu.ucsc.neurobiology.vision.io.PhysiologicalImagingFile(pathToEi);
ei = eiFile.getImage(neuronID); %gets ei data for neuron, storing the information as a 3D array: 2 x nElectrodes x nSamples
clear eiFile

eiData = squeeze(ei(1,2:end,:));

figure
plot(eiData(55,10:40))

plot(eiData(47,10:40))

plot(eiData(26,10:40))

plot(eiData(23,10:40))