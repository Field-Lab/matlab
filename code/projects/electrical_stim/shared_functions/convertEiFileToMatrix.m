function [eiM,neuronIdList] = convertEiFileToMatrix(eiFilePath)
% CONVERTEIFILETOMATRIX(...) takes the path to a .ei file and returns the
% list of neurons and a matrix containing EI waveforms for all neurons
%        input:     eiFilePath is a string ('/path/to/eiFile.ei')
%       ouputs:     eiM: matrix containing EI waveofmrs
%                   neuronIdList: list of neuron IDs
% L Grosberg 9/2015

ei = edu.ucsc.neurobiology.vision.io.PhysiologicalImagingFile(eiFilePath);
neuronIdList = ei.getIDList; 
% Load 1 EI first to allocate memory
neuronEI = ei.getImage(neuronIdList(1));
eiM = zeros(size(neuronEI,2)-1,size(neuronEI,3),size(neuronIdList,1)); 
eiM(:,:,1) = squeeze(neuronEI(1,2:end,:));
 
for n = 2:length(neuronIdList)
    neuronEI = ei.getImage(neuronIdList(n));
    eiM(:,:,n) = squeeze(neuronEI(1,2:end,:)); % Actual EI
end
end