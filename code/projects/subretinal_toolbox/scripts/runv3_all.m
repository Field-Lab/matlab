clear;

if isunix
    visionPath = '/home/ggoetz/Research/Vision/Vision815/Vision.jar'; % Path to Vision - Unix 
    visionWritePath = '/home/ggoetz/Research/Eclipse/110314 - Write Data V4/WriteDataFile.jar';
else
    visionPath = '\\badger\Users\ggoetz\Research\Vision\Vision815\Vision.jar'; 
    visionWritePath = '\\badger\Users\ggoetz\Research\Eclipse\110314 - Write Data V4\WriteDataFile.jar'; 
end

javaaddpath(visionWritePath);
javaaddpath(visionPath);


%% Artifact removal

% Path to the folder in which the data000,... subfolders can be found
rawDataRootPath = '/media/MEA_RAW_6/2013-09-11-0/data/';
% path to the folder in which the data000.txt,... logfiles can be found
logfileRootPath = '/media/MEA_PROCESSED_6/2013-09-11-0/logfiles/';
% path to the output folder. data000,... subfolders will be created there.
processedDataRootPath = '/media/MEA_PROCESSED_6/2013-09-11-0/data/';
 
% Data files that you want to process
dataStrs = {'001','002','003','004','005'};

for kk=1:length(dataStrs)
    rawDataPath = [rawDataRootPath 'data' dataStrs{kk}];
    logfilePath = [logfileRootPath 'logfile' dataStrs{kk} '.txt'];
    processedDataPath = [processedDataRootPath 'data' dataStrs{kk}];
    
    % required parameters
    processData(rawDataPath, logfilePath, processedDataPath,...
                ... % optional parameters under here - see processData help
                            'nTrials',100,... % number of trials used to estimate artifact
                            'nSkipProcessing',20,... % numbers of trials skipped when removing artifact
                            'nPSkipEstimation',25,... % number of trials skipped to estimate artifact
                            'saveArtifact',true,... % save artifacts
                            'visionPath',visionWritePath,... % path to the jar file that has RawDataFile and ModifyRawDataFile classes
                            'smootheData',true); % blanks out part of the artifact
end

% Path to the folder in which the data000,... subfolders can be found
rawDataRootPath = '/media/MEA_RAW_6/2013-09-12-0/data/';
% path to the folder in which the data000.txt,... logfiles can be found
logfileRootPath = '/media/MEA_PROCESSED_6/2013-09-12-0/logfiles/';
% path to the output folder. data000,... subfolders will be created there.
processedDataRootPath = '/media/MEA_PROCESSED_6/2013-09-12-0/data/';
 
% Data files that you want to process
dataStrs = {'001','002','003','004','005','006'};

for kk=1:length(dataStrs)
    rawDataPath = [rawDataRootPath 'data' dataStrs{kk}];
    logfilePath = [logfileRootPath 'logfile' dataStrs{kk} '.txt'];
    processedDataPath = [processedDataRootPath 'data' dataStrs{kk}];
    
    % required parameters
    processData(rawDataPath, logfilePath, processedDataPath,...
                ... % optional parameters under here - see processData help
                            'nTrials',100,... % number of trials used to estimate artifact
                            'nSkipProcessing',20,... % numbers of trials skipped when removing artifact
                            'nPSkipEstimation',25,... % number of trials skipped to estimate artifact
                            'saveArtifact',true,... % save artifacts
                            'visionPath',visionWritePath,... % path to the jar file that has RawDataFile and ModifyRawDataFile classes
                            'smootheData',true); % blanks out part of the artifact
end


%% Neuron statistics computation and plotting

% % % path to the folder in which the data000.txt,... logfiles can be found
% logfileRootPath = '/media/MEA_PROCESSED_6/2013-07-25-0/logfiles/';
% % path to the output folder. data000,... subfolders will be created there.
% processedDataRootPath = '/media/MEA_PROCESSED_6/2013-07-25-0/data/';
% 
% % Data files that you want to process
% dataStrs = {'001','002'};
% 
% for kk=1:length(dataStrs)
%     logfilePath = [logfileRootPath 'logfile' dataStrs{kk} '.txt'];
%     processedDataPath = [processedDataRootPath 'data' dataStrs{kk}];
%     
%     % Neuron statistics
%     analyzeSpikes(processedDataPath,logfilePath,...
%                     'binSize',100,...
%                     'outputPath',[processedDataPath '/statistics/all'],...
%                     'neuronDataFolder',[processedDataPath '/vision_processing/data' dataStrs{kk}]);
% 
%     % Plotting
%     createFigures([processedDataPath filesep 'statistics/all'],[processedDataPath filesep 'figures/all'],... );
%                     'logfile',logfilePath,'combineStim',true,'imageFormat','epsc');                
% end
 
%% Activation data computation
% Set the beginning and end of integration time in computeActivationData.m
% if required
 
% statsRootFolder = '/media/MEA_PROCESSED_5/2013-05-22-0/data/';
% logfileRootPath = '/media/MEA_PROCESSED_5/2013-05-22-0/logfiles/';
% resultRootFolder = '/media/MEA_PROCESSED_5/2013-05-22-0/data/';
% 
% % Setting processing parameters
% powerScalingFactor = 1;       % Used to normalize the powers in the logfile
% powerToIrradianceFactor = 14.68;     % Multiply by max irradiance 
% 
% % Data files that you want to process
% dataStrs = {'001','002'};
% 
%  
% for kk=1:length(dataStrs)
%     logfilePath = [logfileRootPath 'logfile' dataStrs{kk} '.txt'];
%     statsDataPath = [statsRootFolder 'data' dataStrs{kk} '/statistics/all'];
%     resultPath = [resultRootFolder 'data' dataStrs{kk} '/statistics/activation'];
%     
%     computeActivationData(statsDataPath,logfilePath,resultPath,...
%         powerScalingFactor,powerToIrradianceFactor,'createtextfiles',true);
% end

