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
rawDataRootPath = '/media/MEA_RAW_8/2014-01-17-0/data/';
% path to the folder in which the data000.txt,... logfiles can be found
logfileRootPath = '/media/MEA_PROCESSED_8/2014-01-17-0/logfiles/';
% path to the output folder. data000,... subfolders will be created there.
processedDataRootPath = '/media/MEA_PROCESSED_8/2014-01-17-0/data/';
 
% Data files that you want to process
dataStrs = {'001','002','003','004','005'};

for kk=1:length(dataStrs)
    rawDataPath = [rawDataRootPath 'data' dataStrs{kk}];
    logfilePath = [logfileRootPath 'logfile' dataStrs{kk} '.txt'];
    processedDataPath = [processedDataRootPath 'data' dataStrs{kk}];
    
    % required parameters
    processGratingsData(rawDataPath, logfilePath, processedDataPath,...
                ... % optional parameters under here - see processData help
                            'nTrials', 70,... % number of trials used to estimate artifact
                            'ttlsPerPulse', 2,...
                            'nSkipProcessing', 5,... % numbers of trials skipped when removing artifact (NOT TTL PULSES)
                            'nPSkipEstimation', 5,... % number of trials skipped to estimate artifact (NOT TTL PULSES)
                            'saveArtifact', true,... % save artifacts
                            'visionPath', visionWritePath,... % path to the jar file that has RawDataFile and ModifyRawDataFile classes
                            'smootheData', true,... % blanks out part of the artifact
                            'noiseSigma', 5);  % Variance of the WN used when blanking
end

