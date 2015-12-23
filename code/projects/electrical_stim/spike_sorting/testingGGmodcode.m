%%
if ~exist('edu/ucsc/neurobiology/vision/io/RawDataFile','class')
    javaaddpath(visionPath);
end
%%
system = 'stim512'; %'stim64'
pathToLabViewOutputFiles = '/Volumes/Data/2014-08-13-0/'; 
newFileFolder = '/Volumes/Analysis/2014-08-13-0/data003_c';
if ~exist(newFileFolder,'dir')
    mkdir(newFileFolder);
end

datarun = '003';
originalDataFilePath= ['/Volumes/Data/2014-08-13-0/data' datarun '_copy'];
originalRawFile = edu.ucsc.neurobiology.vision.io.RawDataFile(originalDataFilePath);
refHeader = originalRawFile.getHeader();
totalSamples = refHeader.getNumberOfSamples;
Channels = 1:(refHeader.getNumberOfElectrodes-1); %-1 because first channel has TTL pulses
NS_GlobalConstants = NS_GenerateGlobalConstants(length(Channels));

newFile = edu.ucsc.neurobiology.vision.io.ModifyRawDataFile(newFileFolder,refHeader); 


%% 
repeatPeriod = 10000; 

% Copying data to test. 
nSamplesToCopy = totalSamples;
startSample = 0;
while nSamplesToCopy>0
    rawData = originalRawFile.getData(startSample,min(10000,nSamplesToCopy));
    newFile.appendDataToLastFile(rawData);
    nSamplesToCopy = nSamplesToCopy - 10000;
    startSample = startSample + 10000;
end

newFile.close; 
originalRawFile.close; 
clear nSamplesToCopy startSample
%% test files

originalRawFile = edu.ucsc.neurobiology.vision.io.RawDataFile(originalDataFilePath);
newRawFile = edu.ucsc.neurobiology.vision.io.RawDataFile(newFileFolder); 

originalRawData = originalRawFile.getData(0,10000);
newRawData = newRawFile.getData(0,10000); 
if all(originalRawData(:) == newRawData(:))
    disp('Data identical');
end
