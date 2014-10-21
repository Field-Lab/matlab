function allDifferences = compareSingleSpotsToDoubleSpot(singleSpotStatsFolder1, ...
    singleSpotRun1, singleSpotStatsFolder2, singleSpotRun2, ...
    doubleSpotStatsFolder, doubleSpotRun)
% This function computes the difference in number of spikes between
% different stimulation runs. It is meant to investigate possible
% non-linearity effects of double spot stimulation.
%
% Parameters:
%
% 
% Returns:
%
%
% Version: 0.1 - 2013/04/24
%

nBins = 50;

%% Load activation data for spot 1, spot 2 and the double spot

% Finding the neuron IDs for each dataset, final list is those in common
neuronIDsSpot1 = findNeuronIDsInStatsFolder(singleSpotStatsFolder1);
neuronIDsSpot2 = findNeuronIDsInStatsFolder(singleSpotStatsFolder2);
neuronIDsDoubleSpot = findNeuronIDsInStatsFolder(doubleSpotStatsFolder);

neuronIDs = intersect(neuronIDsSpot1, intersect(neuronIDsSpot2, neuronIDsDoubleSpot));
nNeurons = length(neuronIDs);

% Loading the activation data of each of those neurons for the relevant
% stimulus run
actData = zeros(nNeurons,3);
load(fullfile(singleSpotStatsFolder1, sprintf('neuron%d_act.mat',neuronIDs(1))));
stimDataPos = ismember(activationData.labels,'nSpikes');

for kk=1:nNeurons
    % Spot 1...
    load(fullfile(singleSpotStatsFolder1, sprintf('neuron%d_act.mat',neuronIDs(kk))));
    actData(kk,1) = activationData.data(singleSpotRun1,stimDataPos);
    % Spot 2...
    load(fullfile(singleSpotStatsFolder2, sprintf('neuron%d_act.mat',neuronIDs(kk))));
    actData(kk,2) = activationData.data(singleSpotRun2,stimDataPos);
    % Double spot
    load(fullfile(doubleSpotStatsFolder, sprintf('neuron%d_act.mat',neuronIDs(kk))));
    actData(kk,3) = activationData.data(doubleSpotRun,stimDataPos);
end


%% Compare this activation data

allDifferences = actData(:,3) - actData(:,1) - actData(:,2);


%% Create plots

fh = figure(); set(fh,'color','white');
hist(allDifferences, nBins);
title(sprintf('%s; spot1 = %d, spot2 = %d, double spot = %d',...
    'Spikes elicited by double spot vs sum of single spot',...
    singleSpotRun1, singleSpotRun2, doubleSpotRun));
xlabel('Additional number of spikes in double spot')
ylabel('Number of neurons')

fh = figure(fh+1); set(fh,'color','white');
hold on
scatter(actData(:,1)+actData(:,2),actData(:,3));
plot(-2:5:15,-2:5:15,'r');
title(sprintf('%s; spot1 = %d, spot2 = %d, double spot = %d',...
    'Spikes elicited by double spot vs sum of single spot',...
    singleSpotRun1, singleSpotRun2, doubleSpotRun));
xlabel('Sum of spikes in single spot')
ylabel('Number of spikes in double spot')
axis([0 8 0 8])

end % compateSingleSpotsToDoubleSpot


function neuronIDList = findNeuronIDsInStatsFolder(statsFolder)
% This function returns a list of all neuron IDs in a given stats folder.

% Finding the neuron IDs, spot 1
contentsStatsFolder = dir(statsFolder);
neuronIDList = [];
for kk=1:length(contentsStatsFolder)
    pos = strfind(contentsStatsFolder(kk).name,'.mat');
    if (pos)
        neuronIDList(end+1) = str2double(contentsStatsFolder(kk).name(7:(pos-5)));
    end
end

end % findNeuronIDsInStatsFolder