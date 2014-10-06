function [amps channelsWithStim stimAmpVectors channelsConnected] = getStimAmps(dataPath, patternNumber, movieNumber)

% arguments
%   dataPath: path from current directory to pattern definitions file
%
%
%

if isnumeric(patternNumber)
    patternNumber = num2str(patternNumber);
end

NS_GlobalConstants=NS_GenerateGlobalConstants(61);

if exist([dataPath filesep 'status_files'], 'file') 
    FullName_status = [dataPath filesep 'status_files' filesep 'status_m' num2str(movieNumber)];
else
    FullName_status=[dataPath filesep 'status_m' num2str(movieNumber)];
end

if ~exist([FullName_status '.mat'], 'file')
   warnH = warndlg([FullName_status '.mat could not be found']);
   uiwait(warnH)
   return
end
load(FullName_status); %contains Status
if ~exist('status', 'var') %because variable status is only sometimes capitalized
    status = Status;
end


if exist([dataPath filesep 'pattern_files'], 'file') 
    FullName_pattern = [dataPath filesep 'pattern_files' filesep 'pattern' patternNumber '_m' num2str(movieNumber)];
else
    FullName_pattern = [dataPath filesep 'pattern' patternNumber '_m' num2str(movieNumber)];
end

if ~exist([FullName_pattern '.mat'], 'file')
   warnH = warndlg([FullName_pattern '.mat could not be found']);
   uiwait(warnH)
   return
end

load(FullName_pattern); %contains Pattern
%size(Pattern)
if ~exist('Pattern', 'var') %because variable pattern is only sometimes capitalized
    Pattern = pattern;
end
    
thisPatternLong = Pattern;

% determines whether any secondary electrodes were "connected" (even though the pulse amplitude may have rounded to zero)
channelsConnected = [];
for ii = 1:length(thisPatternLong)
    if any(thisPatternLong(ii).data(3,:) ~= 0); %if "connect" signal is nonzero at some point
        channelsConnected = [channelsConnected thisPatternLong(ii).channel];
    end
end

%remove elements of thisPattern in which no stimulation occurs
thisPattern = thisPatternLong;
for i = length(thisPatternLong):-1:1 %looped backwards so that removal of elements works correctly
    if min(thisPatternLong(i).data(1,:) == 0) == 1; %if all DAC values for this electrode = 0
        thisPattern(i) = [];
    end
end

nElecInPattern = length(thisPattern);

%get the current ranges for each electrode included in the pattern
currentRanges = NS_GlobalConstants.CurrentRanges;

elecCurrentStep = zeros(1, nElecInPattern);
for i = 1:nElecInPattern
    elecCurrentStep(i) = currentRanges(status.ChannelsStatus(thisPattern(i).channel).range + 1)/127;
end

%get channel numbers from thisPattern and finds where they are within argument "channels"
channelsWithStim = zeros(1, nElecInPattern);
for i = 1:nElecInPattern
    channelsWithStim(i) = thisPattern(i).channel;
end

%create vector of different current amplitudes used in pattern
patternLength = size(thisPattern(1).data,2);
stimAmpVectors = zeros(nElecInPattern, patternLength);
amps = zeros(nElecInPattern, 1);

for i = 1:nElecInPattern
    for j = 1:patternLength
        stimAmpVectors(i,j) = thisPattern(i).data(1,j)*thisPattern(i).data(3,j)*elecCurrentStep(i);
    end
    stimAmpAbsMax = max(max(abs(stimAmpVectors(i,:))));
    maxIndex = find(stimAmpAbsMax == squeeze(abs(stimAmpVectors(i,:))), 1);
    amps(i) = stimAmpVectors(i, maxIndex);
end
