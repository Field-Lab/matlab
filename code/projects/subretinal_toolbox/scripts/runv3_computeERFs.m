clear;

deviceType = 'small';
neuronList1 = [196,212,227,286,301,691,...
    5206,5567,5673,6121,6376,6557,6601,...
    6751,6978,7127,7156,7621]; % Neurons stimulated with FF
neuronList2 = []; % Neurons stimulated with single px
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

% pixelOrder = [1,3,2,...
%               9,6,5,4,7,8,...
%               15,12,10,14,13,16,18,11,17,...
%               21,20,24,25,23,28,27,19,22,29,26,...
%               37,31,41,36,38,40,33,32,35,30,39,34,...
%               48,49,43,51,50,47,45,46,42,53,52,44,...
%               64,62,60,57,63,56,61,58,65,59,54,55,...
%               68,72,76,69,74,70,75,71,67,66,73,77,...
%               81,86,87,79,89,83,84,88,80,82,85,78,...
%               94,95,90,96,99,93,101,100,97,91,92,98,...
%               104,102,107,112,111,110,105,109,103,108,106,113,...
%               120,114,123,116,118,115,124,119,121,117,122,...
%               126,133,129,130,127,128,125,132,131,...
%               137,138,139,134,135,136,...
%               141,142,140];
          
pixelOrder = [3,1,2,6,7,8,5,4,9,16,11,17,14,15,18,13,10,12,22,28,29,21,...
    25,27,26,19,23,24,20,32,40,34,35,30,33,38,41,37,31,36,39,45,51,42,47,...
    46,48,53,52,49,50,43,44,55,54,60,62,58,56,63,59,64,61,65,57,68,69,77,...
    74,72,70,76,67,73,71,75,66,78,86,81,89,80,84,83,79,88,87,82,85,99,92,...
    91,96,101,100,98,93,90,97,94,95,112,105,113,107,103,111,104,106,108,...
    109,110,102,123,119,121,115,120,117,124,116,114,122,118,127,133,129,...
    128,131,130,125,126,132,136,135,134,139,138,137,140,142,141];

% columns 1-5
pixelIndex = 1:41;
psthStatsFolder1 = '/media/MEA_PROCESSED_7/2013-10-09-0/data/data001/statistics/all';
activationStatsFolder1 = '/media/MEA_PROCESSED_7/2013-10-09-0/data/data001/statistics/activation_nosupp';
outputFolderStatistics1 = '/media/MEA_PROCESSED_7/2013-10-09-0/data/data001/statistics/e_rfs/stim_neurons_nosupp';

computeERFs(psthStatsFolder1, activationStatsFolder1, deviceType, ...
    pixelOrder(pixelIndex), outputFolderStatistics1, ...
    'neuronList',neuronList); ...,'neuronPosition',neuronPos);

% columns 6-10
pixelIndex = 42:101;
psthStatsFolder2 = '/media/MEA_PROCESSED_7/2013-10-09-0/data/data002/statistics/all';
activationStatsFolder2 = '/media/MEA_PROCESSED_7/2013-10-09-0/data/data002/statistics/activation_nosupp';
outputFolderStatistics2 = '/media/MEA_PROCESSED_7/2013-10-09-0/data/data002/statistics/e_rfs/stim_neurons_nosupp';

computeERFs(psthStatsFolder2, activationStatsFolder2, deviceType, ...
    pixelOrder(pixelIndex), outputFolderStatistics2, ...
    'neuronList',neuronList); ...,'neuronPosition',neuronPos);

% columns 11-15
pixelIndex = 102:142;
psthStatsFolder3 = '/media/MEA_PROCESSED_7/2013-10-09-0/data/data003/statistics/all';
activationStatsFolder3 = '/media/MEA_PROCESSED_7/2013-10-09-0/data/data003/statistics/activation_nosupp';
outputFolderStatistics3 = '/media/MEA_PROCESSED_7/2013-10-09-0/data/data003/statistics/e_rfs/stim_neurons_nosupp';

computeERFs(psthStatsFolder3, activationStatsFolder3, deviceType, ...
    pixelOrder(pixelIndex), outputFolderStatistics3, ...
    'neuronList',neuronList); ...,'neuronPosition',neuronPos);



%% Combine all the ERFs

outputFolderStatistics = '/media/MEA_PROCESSED_7/2013-10-09-0/results/erfs/statistics/stim_neurons_nosupp';
allSourceFolders = {outputFolderStatistics1,outputFolderStatistics2,outputFolderStatistics3};

combineERFs(outputFolderStatistics,allSourceFolders);

%% Plot the eRFs

outputFolderFigures = '/media/MEA_PROCESSED_7/2013-10-09-0/results/delete_me';
plotERFs(outputFolderStatistics, outputFolderFigures);

% outputFolderFigures = '/media/MEA_PROCESSED_7/2013-10-09-0/results/erfs_nosupp/figures/td';
% plotERFMovie(outputFolderStatistics, outputFolderFigures, 'nFrames', 40);

%% Get the distribution of diameters


nNeurons = length(neuronList);
allDiamsEL = zeros(nNeurons,1);
for kk=1:nNeurons
    load(fullfile(outputFolderStatistics, sprintf('neuron%d_erf.mat',neuronList(kk))));
    eRF = 1./sqrt(svd(neuronERFdata.sigma));
    allDiamsEL(kk) = 2*geo_mean(eRF);
end

