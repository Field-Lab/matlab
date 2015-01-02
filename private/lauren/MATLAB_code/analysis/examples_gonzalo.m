%% Get electrode positions. 

[xCoords, yCoords] = getElectrodeCoords512(); 
positions = loadElecPositions512();
figure; scatter(xCoords, yCoords,1,'white'); axis image; axis off; 
for e = 1:512
    text(positions(e,1),positions(e,2),num2str(e),'HorizontalAlignment','center');
end



%% Load electrophysiological image (EI) to get cell templates

neuronPath = '/Volumes/Analysis/2014-11-05-8/data001/data001.neurons';
neuronFile = edu.ucsc.neurobiology.vision.io.NeuronFile(neuronPath);
neuronIds = neuronFile.getIDList(); %returns a list of all the neuron IDs in this neuron file

ei_path = '/Volumes/Analysis/2014-11-05-8/data001/data001.ei';
eiFile = edu.ucsc.neurobiology.vision.io.PhysiologicalImagingFile(ei_path);
% returns a matrix representing the electrophysiological image of the 
% neuron with ID neuronId. This matrix has a dimension of <math>2 \times 
% nElectrodes \times depth</math>, where <math>nElectrodes</math> is the 
% number of electrodes of the array and <math>depth</math> the depth of 
% the EI. Then, in Matlab syntax, ei(1,:,:) corresponds to the value of 
% the EI, and ei(2,:,:) corresponds to the variance in the EI.

neuronId = neuronIds(101); 
neuronEI = eiFile.getImage(neuronId); 
neuronEI_volt = squeeze(neuronEI(1,:,:)).'; % Actual EI
neuronEI_var = squeeze(neuronEI(2,:,:)).';  % Variances

% Plot the EI
figure; set(gcf,'Color','white'); 
arbitraryScaling = 10; 
scatter(xCoords, yCoords, abs(min(neuronEI_volt(:,2:end)))*arbitraryScaling);
axis image; axis off; title(sprintf('EI for neuron %0.0f',neuronId)); 

overlayElecNums = 1; 
if overlayElecNums
    for e = 1:512
        text(positions(e,1),positions(e,2),num2str(e),'HorizontalAlignment','center');
    end
end
%% Look at the waveforms on a cluster of electrodes

% Get cluster of electrodes
centerElectrode   = 167; %example; 
clusterElectrodes = getCluster512(centerElectrode);

figure; plot(neuronEI_volt(:,clusterElectrodes)); 
xlabel('samples'); ylabel('amplitude'); 
title('templates on a cluster of electrodes'); 

