clear;

experimentDate = '2013-10-16-0';
deviceType = 'small';
neuronList1 = [2659,2881,2896,6106,6346,...
    6511,6721]; % Neurons stimulated with FF
neuronList2 = [2297,2581,2642,2716,2896,3272,6106,6346,6496,6857]; % Neurons stimulated with single px
neuronList = union(neuronList1, neuronList2);

%% Map the coordinates
% We want to know where the pixels are on the MEA and where the neurons are
% on the pixels

% pixelIndex = [1,19,41,90,113,142];
% pixelCoordsOnMEA = [-225   430;
%                      105   450;
%                     -495   -70;
%                      390   180;
%                     -255  -390;
%                      180  -360];
% flipped = false;
% dataFolder = '/media/MEA_PROCESSED_7/2013-10-09-0/data/data001/vision_processing/data001';
% 
% [arrayToMEA MEAToArray] = mapCoordinates(pixelCoordsOnMEA, pixelIndex, deviceType, flipped);
% [neuronPos neuronList] = computeNeuronsCoordinatesOnPVArray(dataFolder, MEAToArray, neuronList);

%% Compute and plot the electric receptive fields

pixelOrder = [2,3,1,6,8,7,4,9,5,10,14,16,15,18,11,17,12,13,23,27,22,20,...
    21,26,25,29,28,24,19,39,33,35,36,34,37,30,32,38,31,40,41,52,47,53,50,...
    44,46,49,51,43,42,48,45,61,62,57,63,58,59,56,64,65,60,54,55,77,73,71,...
    69,74,67,68,75,66,72,70,76,83,89,81,82,80,84,85,88,79,87,86,78,98,101,...
    99,92,100,91,95,97,96,90,94,93,105,111,103,104,106,108,110,107,102,...
    112,109,113,115,123,122,116,119,124,117,114,118,120,121,131,128,129,...
    125,126,127,133,130,132,136,137,134,139,135,138,141,142,140];

% columns 1-5
pixelIndex = 1:41;
psthStatsFolder1 = sprintf('/media/MEA_PROCESSED_7/%s/data/data001/statistics/all', experimentDate);
activationStatsFolder1 = sprintf('/media/MEA_PROCESSED_7/%s/data/data001/statistics/activation_nosupp', experimentDate);
outputFolderStatistics1 = sprintf('/media/MEA_PROCESSED_7/%s/data/data001/statistics/e_rfs/stim_neurons', experimentDate);

computeERFs(psthStatsFolder1, activationStatsFolder1, deviceType, ...
    pixelOrder(pixelIndex), outputFolderStatistics1, ...
    'neuronList',neuronList); ...,'neuronPosition',neuronPos);

% columns 6-10
pixelIndex = 42:101;
psthStatsFolder2 = sprintf('/media/MEA_PROCESSED_7/%s/data/data002/statistics/all', experimentDate);
activationStatsFolder2 = sprintf('/media/MEA_PROCESSED_7/%s/data/data002/statistics/activation_nosupp', experimentDate);
outputFolderStatistics2 = sprintf('/media/MEA_PROCESSED_7/%s/data/data002/statistics/e_rfs/stim_neurons', experimentDate);

computeERFs(psthStatsFolder2, activationStatsFolder2, deviceType, ...
    pixelOrder(pixelIndex), outputFolderStatistics2, ...
    'neuronList',neuronList); ...,'neuronPosition',neuronPos);

% columns 11-15
pixelIndex = 102:142;
psthStatsFolder3 = sprintf('/media/MEA_PROCESSED_7/%s/data/data003/statistics/all', experimentDate);
activationStatsFolder3 = sprintf('/media/MEA_PROCESSED_7/%s/data/data003/statistics/activation_nosupp', experimentDate);
outputFolderStatistics3 = sprintf('/media/MEA_PROCESSED_7/%s/data/data003/statistics/e_rfs/stim_neurons', experimentDate);

computeERFs(psthStatsFolder3, activationStatsFolder3, deviceType, ...
    pixelOrder(pixelIndex), outputFolderStatistics3, ...
    'neuronList',neuronList); ...,'neuronPosition',neuronPos);



%% Combine all the ERFs

outputFolderStatistics = sprintf('/media/MEA_PROCESSED_7/%s/results/erfs/statistics_nosupp/stim_neurons', experimentDate);
allSourceFolders = {outputFolderStatistics1,outputFolderStatistics2,outputFolderStatistics3};

combineERFs(outputFolderStatistics,allSourceFolders);

%% Plot the eRFs

sprintf('/media/MEA_PROCESSED_7/%s/results/erfs_nosupp/figures/still/eps', experimentDate);
outputFolderFigures = sprintf('/media/MEA_PROCESSED_7/%s/results/erfs/figures_nosupp/still/eps', experimentDate);
plotERFs(outputFolderStatistics, outputFolderFigures);

% outputFolderFigures = sprintf('/media/MEA_PROCESSED_7/%s/results/erfs_nosupp/figures/td', experimentDate);
% plotERFMovie(outputFolderStatistics, outputFolderFigures, 'nFrames', 40);

%% Get the distribution of diameters


% nNeurons = length(neuronList);
% allDiamsEL = zeros(nNeurons,1);
% for kk=1:nNeurons
%     load(fullfile(outputFolderStatistics, sprintf('neuron%d_erf.mat',neuronList(kk))));
%     eRF = 1./sqrt(svd(neuronERFdata.sigma));
%     allDiamsEL(kk) = 2*geo_mean(eRF);
% end

