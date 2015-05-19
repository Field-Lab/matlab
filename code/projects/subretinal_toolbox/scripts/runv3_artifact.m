clear;
dbstop if error;

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
rawDataRootPath = '/media/MEA_RAW_8/2014-03-03-0/data/';
% path to the folder in which the data000.txt,... logfiles can be found
logfileRootPath = '/media/MEA_PROCESSED_8/2014-03-03-0/logfiles/';
% path to the output folder. data000,... subfolders will be created there.
processedDataRootPath = '/media/MEA_PROCESSED_8/2014-03-03-0/data/';
 
% Data files that you want to process
dataStrs = {'000'};

ttlPulseTimes = {(0:150)*200000+20*20000+106};

for kk=1:length(dataStrs)
    rawDataPath = [rawDataRootPath 'data' dataStrs{kk}];
    logfilePath = [logfileRootPath 'logfile' dataStrs{kk} '.txt'];
    processedDataPath = [processedDataRootPath 'data' dataStrs{kk}];
    
    % required parameters
    processData(rawDataPath, logfilePath, processedDataPath,...
                ... % optional parameters under here - see processData help
                            'nTrials',100,... % number of trials used to estimate artifact
                            'nSkipProcessing',5,... % numbers of trials skipped when removing artifact
                            'nPSkipEstimation',10,... % number of trials skipped to estimate artifact
                            'saveArtifact',true,... % save artifacts
                            'visionPath',visionWritePath,... % path to the jar file that has RawDataFile and ModifyRawDataFile classes
                            'smootheData',false,... % blanks out part of the artifact
                            'noiseSigma',0.2,... % Variance of the WN used when blanking
                            'logfileType','labview',... % Logfile generation mechanism
                            'ttlTimes',ttlPulseTimes);  
end

