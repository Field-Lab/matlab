clear;

deviceType = 'small';
neuronList1 = []; % Neurons stimulated with FF
neuronList2 = [1006 1112 1156 1246 1306 1352 1486 1546 1651 1787 31 ...
    691 706 7561 76 796 871 931]; % Neurons stimulated with single px
neuronList = union(neuronList1, neuronList2);

%% Map the coordinates
% We want to know where the pixels are on the MEA and where the neurons are
% on the pixels

pixelIndex = [19,90,113,142];
pixelCoordsOnMEA = [ 465 -390;
                     315   30; 
                    -405 -370; 
                    -285   30];
flipped = false;
dataFolder = '/media/MEA_PROCESSED_6/2013-07-09-0/data/data002/vision_processing/data002';

[arrayToMEA MEAToArray] = mapCoordinates(pixelCoordsOnMEA, pixelIndex, deviceType, flipped);
[neuronPos neuronList] = computeNeuronsCoordinatesOnPVArray(dataFolder, MEAToArray, neuronList);

%% Compute and plot the electric receptive fields

% columns 6-10
pixelIndex = 42:101;
activationStatsFolder1 = '/media/MEA_PROCESSED_6/2013-07-09-0/data/data002/statistics/activation';
outputFolderStatistics1 = '/media/MEA_PROCESSED_6/2013-07-09-0/data/data002/statistics/e_rfs/stimulated_neurons_with_eis';
outputFolderFigures1 = '/media/MEA_PROCESSED_6/2013-07-09-0/data/data002/figures/e_rfs/stimulated_neurons_with_eis';

computeElectricReceptiveField(activationStatsFolder1, deviceType, ...
    pixelIndex, outputFolderStatistics1, outputFolderFigures1, ...
    'neuronList',neuronList,'savePlots',false,'neuronPosition',neuronPos);

% columns 11-15
pixelIndex = 102:142;
activationStatsFolder2 = '/media/MEA_PROCESSED_6/2013-07-09-0/data/data003/statistics/activation';
outputFolderStatistics2 = '/media/MEA_PROCESSED_6/2013-07-09-0/data/data003/statistics/e_rfs/stimulated_neurons_with_eis';
outputFolderFigures2 = '/media/MEA_PROCESSED_6/2013-07-09-0/data/data003/figures/e_rfs/stimulated_neurons_with_eis';

computeElectricReceptiveField(activationStatsFolder2, deviceType, ...
    pixelIndex, outputFolderStatistics2, outputFolderFigures2, ...
    'neuronList',neuronList,'savePlots',false,'neuronPosition',neuronPos);

% columns 1-5
pixelIndex = 1:41;
activationStatsFolder3 = '/media/MEA_PROCESSED_6/2013-07-09-0/data/data004/statistics/activation';
outputFolderStatistics3 = '/media/MEA_PROCESSED_6/2013-07-09-0/data/data004/statistics/e_rfs/stimulated_neurons_with_eis';
outputFolderFigures3 = '/media/MEA_PROCESSED_6/2013-07-09-0/data/data004/figures/e_rfs/stimulated_neurons_with_eis';

computeElectricReceptiveField(activationStatsFolder3, deviceType, ...
    pixelIndex, outputFolderStatistics3, outputFolderFigures3, ...
    'neuronList',neuronList,'savePlots',false,'neuronPosition',neuronPos);


%% Combine all the ERFs

outputFolderFigures = '/media/MEA_PROCESSED_6/2013-07-09-0/results/erfs/figures_erfs_combined_with_eis';
outputFolderStatistics = '/media/MEA_PROCESSED_6/2013-07-09-0/results/erfs/statistics_combined_with_eis';
allSourceFolders = {outputFolderStatistics1,outputFolderStatistics2,outputFolderStatistics3};

combineAndPlotERFs(outputFolderFigures,outputFolderStatistics,allSourceFolders,neuronPos);

%% Get the distribution of diameters

% Col1: ERF ID; Col2: WN ID
idCorrespondance = [31   91;
                    76   76;
                    691  691;
                    706  706;
                    796  796;
                    871  871;
                    931  931;
                    1006 1006;
                    1112 1111;
                    1156 1156;
                    1246 1246;
                    1306 1306;
                    1352 1231;
                    1486 1486;
                    1546 1546;
                    1651 1651;
                    1787 -1;
                    7561 167];

nNeurons = size(idCorrespondance,1);
allDiamsEL = zeros(nNeurons,1);
for kk=1:nNeurons
    load(fullfile(outputFolderStatistics, sprintf('neuron%d_erf.mat',idCorrespondance(kk,1))));
    eRF = 1./sqrt(svd(neuronERFdata.sigma));
    allDiamsEL(kk) = 2*geo_mean(eRF);
end

% Get the same information for the visible light data
% Compute the diameters for the WN dataset

paramsFileWNPath = '/media/MEA_PROCESSED_6/2013-07-09-0/data/data000/vision_processing/data000/data000.params';
paramsFileWN = edu.ucsc.neurobiology.vision.io.ParametersFile(paramsFileWNPath);
allDiamsWN = zeros(nNeurons,1);
for kk=1:nNeurons
    cID = idCorrespondance(kk,2);
    if cID>0
        cRF(1) = paramsFileWN.getDoubleCell(cID,'SigmaX');
        cRF(2) = paramsFileWN.getDoubleCell(cID,'SigmaY');
        allDiamsWN(kk) = 2*geo_mean(cRF);
    else
        allDiamsWN(kk) = NaN;
    end
end