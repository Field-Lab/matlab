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


%% Neuron statistics computation and plotting

% % path to the folder in which the data000.txt,... logfiles can be found
logfileRootPath = '/media/MEA_PROCESSED_7/2013-12-03-0/logfiles/';
% path to the output folder. data000,... subfolders will be created there.
processedDataRootPath = '/media/MEA_PROCESSED_7/2013-12-03-0/data/';

% Data files that you want to process
dataStrs = {'003'};

for kk=1:length(dataStrs)
    logfilePath = [logfileRootPath 'logfile' dataStrs{kk} '.txt'];
    processedDataPath = [processedDataRootPath 'data' dataStrs{kk}];
    
    % Neuron statistics
    analyzeSpikes(processedDataPath,logfilePath,...
                    'binSize',100,...
                    'ttlsPerPulse',3,...
                    'outputPath',[processedDataPath '/statistics/all'],...
                    'neuronDataFolder',[processedDataPath '/vision_processing/data' dataStrs{kk}]);

    % Plotting
    createFigures([processedDataPath filesep 'statistics/all'],[processedDataPath filesep 'figures/all'],... );
                    'logfile',logfilePath,'combineStim',true,'imageFormat','epsc');                
end
 
% %% Activation data computation
% % Set the beginning and end of integration time in computeActivationData.m
% % if required
%  
% statsRootFolder = '/media/MEA_PROCESSED_7/2013-10-22-0/data/';
% logfileRootPath = '/media/MEA_PROCESSED_7/2013-10-22-0/logfiles/';
% resultRootFolder = '/media/MEA_PROCESSED_7/2013-10-22-0/data/';
% 
% % Setting processing parameters
% powerScalingFactor = 1;       % Used to normalize the powers in the logfile
% powerToIrradianceFactor = 9.5;     % Multiply by max irradiance 
% 
% % Data files that you want to process
% dataStrs = {'001','002','003','004'};
% 
%  
% for kk=1:length(dataStrs)
%     logfilePath = [logfileRootPath 'logfile' dataStrs{kk} '.txt'];
%     statsDataPath = [statsRootFolder 'data' dataStrs{kk} '/statistics/all'];
%     resultPath = [resultRootFolder 'data' dataStrs{kk} '/statistics/activation'];
%     
%     computeActivationData(statsDataPath,logfilePath,resultPath,...
%         powerScalingFactor,powerToIrradianceFactor,'createtextfiles',true,...
%         'positiveOnly',false);
% end

