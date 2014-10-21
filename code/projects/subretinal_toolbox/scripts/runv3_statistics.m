clear;

if isunix
    visionPath = '/home/ggoetz/Research/Vision/UCSC/Vision_8.2.4/Vision.jar'; % Path to Vision - Unix 
    visionWritePath = '/home/ggoetz/Research/Eclipse/110314 - Write Data V4/WriteDataFile.jar';
else
    visionPath = '\\badger\Users\ggoetz\Research\Vision\Vision815\Vision.jar'; 
    visionWritePath = '\\badger\Users\ggoetz\Research\Eclipse\110314 - Write Data V4\WriteDataFile.jar'; 
end

javaaddpath(visionWritePath);
javaaddpath(visionPath);

%% Procesing parameters

logfiletype = 'labview';

% path to the output folder. data000,... subfolders will be created there.
processedDataRootPath = '/media/MEA_PROCESSED_9/2014-07-11-0/data/';
% path to the folder in which the data000.txt,... logfiles can be found
logfileRootPath = '/media/MEA_PROCESSED_9/2014-07-11-0/logfiles/alternate_single-pulse';
% all stats folders are subdirectories of this root folder
statsRootFolder = '/media/MEA_PROCESSED_9/2014-07-11-0/data/';
% all figure folders are subdirectories of this root folder
resultRootFolder = '/media/MEA_PROCESSED_9/2014-07-11-0/data/';

% Data files that you want to process
dataStrs = {'001'};

%% Neuron statistics computation and plotting

for kk=1:length(dataStrs)
    logfilePath = fullfile(logfileRootPath, ['logfile' dataStrs{kk} '.txt']);
    processedDataPath = fullfile(processedDataRootPath, ['data' dataStrs{kk}]);
    
    % Neuron statistics
    analyzeSpikes(processedDataPath,logfilePath,...
                    'logfileType',logfiletype,...  
                    'binSize',100,...
                    'ttlsPerPulse',1,...
                    'outputPath',fullfile(processedDataPath, 'statistics/all'),...
                    'neuronDataFolder',fullfile(processedDataPath, ['vision_processing/data' dataStrs{kk}]));

    % Plotting
    createFigures(fullfile(processedDataPath, 'statistics/all'),...
                  fullfile(processedDataPath, 'figures/all'),... 
                  'logfileType',logfiletype,... 
                  'logfile',logfilePath,...
                  'combineStim',true,...
                  'imageFormat','epsc');                
end
 
%% Activation data computation
% Set the beginning and end of integration time in computeActivationData.m
% if required

% Setting processing parameters
powerScalingFactor = 1;             % Used to normalize the powers in the logfile
powerToIrradianceFactor = 15.6;     % Multiply by max irradiance 
 
for kk=1:length(dataStrs)
    logfilePath = fullfile(logfileRootPath,  ['logfile' dataStrs{kk} '.txt']);
    statsDataPath = fullfile(statsRootFolder,  ['data' dataStrs{kk} '/statistics/all']);
    resultPath = fullfile(resultRootFolder,  ['data' dataStrs{kk} '/statistics/activation']);
    
    computeActivationData(statsDataPath,logfilePath,resultPath,...
        powerScalingFactor,powerToIrradianceFactor,...
        'logfileType',logfiletype,...
        'createtextfiles',true,...
        'positiveOnly',false,...
        'useExpIdAsBaseline',1); % Note: experiment IDs are 1-based
end

