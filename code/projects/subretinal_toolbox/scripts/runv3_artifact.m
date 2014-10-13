clear;
dbstop if error;

if isunix
    visionPath = '/home/ggoetz/Research/Vision/UCSC/Vision_8.2.4/Vision.jar'; % Path to Vision - Unix 
    visionWritePath = '/home/ggoetz/Research/Eclipse/110314 - Write Data V4/WriteDataFile.jar';
else
    visionPath = '\\badger\Users\ggoetz\Research\Vision\Vision815\Vision.jar'; 
    visionWritePath = '\\badger\Users\ggoetz\Research\Eclipse\110314 - Write Data V4\WriteDataFile.jar'; 
end

javaaddpath(visionWritePath);
javaaddpath(visionPath);


%% Artifact removal

% Path to the folder in which the data000,... subfolders can be found
rawDataRootPath = '/Volumes/Data/UCSC/2014-10-03-0/data/';
% path to the folder in which the data000.txt,... logfiles can be found
logfileRootPath = '/Volumes/Analysis/UCSC/2014-10-03-0/logfiles/';
% path to the output folder. data000,... subfolders will be created there.
processedDataRootPath = '/Volumes/Analysis/UCSC/2014-10-03-0/data/';

% Data files that you want to process
dataStrs = {'001', '003', '005'};

% % Typical frequency reconstruction
% ttlPulseTimes = [{int32((0:199)*1*20000 + 11*20000 + 106), ...
%                   int32(round((0:198)*1/2*20000) + 216*20000 + 106), ...
%                   int32(round((0:298)*1/3*20000) + 321*20000 + 106), ...
%                   int32(round((0:397)*1/4*20000) + 426*20000 + 106), ...
%                   int32(round((0:497)*1/5*20000) + 531*20000 + 106), ...
%                   int32(round((0:596)*1/6*20000) + 636*20000 + 106), ...
%                   int32(round((0:696)*1/7*20000) + 741*20000 + 106), ...
%                   int32(round((0:795)*1/8*20000) + 846*20000 + 106), ...
%                   int32(round((0:895)*1/9*20000) + 951*20000 + 106), ...
%                   int32(round((0:494)*1/10*20000) + 1056*20000 + 106), ...
%                   int32(round((0:593)*1/12*20000) + 1111*20000 + 106), ...
%                   int32(round((0:692)*1/14*20000) + 1166*20000 + 106), ...
%                   int32(round((0:791)*1/16*20000) + 1221*20000 + 106), ...
%                   int32(round((0:890)*1/18*20000) + 1276*20000 + 106), ...
%                   int32(round((0:989)*1/20*20000) + 1331*20000 + 106)}];

% % Typical contrast reconstruction
ttlPulseTimes = [{(0:148)*10*20000 + 10*20000 + 106}, ...
                 {(0:148)*10*20000 + 10*20000 + 106}, ...
                 {(0:148)*10*20000 + 10*20000 + 106}];

for kk=1:length(dataStrs)
    rawDataPath = fullfile(rawDataRootPath, ['data' dataStrs{kk}]);
    logfilePath = fullfile(logfileRootPath, ['logfile' dataStrs{kk} '.txt']);
    processedDataPath = fullfile(processedDataRootPath, ['data' dataStrs{kk}]);
    
    % required parameters
    processData(rawDataPath, logfilePath, processedDataPath,...
                ... % optional parameters under here - see processData help
                            'nTrials',50,...             % number of trials used to estimate artifact
                            'ttlsPerPulse',1,...         % Number of trigger pulses per artefact
                            'nSkipProcessing',0,...      % numbers of trials skipped when removing artifact
                            'nPSkipEstimation',10,...    % number of trials skipped to estimate artifact
                            'saveArtifact',true,...      % save artifacts
                            'visionPath',visionWritePath,... % path to the jar file that has RawDataFile and ModifyRawDataFile classes
                            'smootheData',true,...       % blanks out part of the artifact
                            'blankOnly',false,...        % do we only want to blank or also subtract things?
                            'noiseSigma',2,...           % Variance of the WN used when blanking
                            'logfileType','labview',...  % Logfile generation mechanism
                            'ttlTimes',ttlPulseTimes(kk),... % Manual specification of TTL pulse times
                            'trigToPulseDelay', 15);     % Delay between TTL and artefact. 2 is a good value unless you *know* you want something else.
end

