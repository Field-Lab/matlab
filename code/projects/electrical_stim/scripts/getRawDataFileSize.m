function [totalSamples totalTime] = getRawDataFileSize(pathToData)

% finds the length of the dataset (defined by folder data***) in samples and time (s)
% argument: 
%    pathToData: path from root to directory containing binary raw data files
%    e.g. /Volumes/Palace/Data/Lauren/2008-08-0-27-4/data000/
% 
% assumes that sampe rate is 20000 Hz
%

%makes sure pathToData ends with / (mac, unix) or \ (pc)
if ~strcmp(pathToData(end), filesep)
    pathToData = [pathToData filesep];
end

rawFile = edu.ucsc.neurobiology.vision.io.RawDataFile(pathToData);

header = rawFile.getHeader();
sampleSize = header.getSampleSize();
headerSize = header.getHeaderSize();

fileNames = dir([pathToData 'data*']);
totalSize = 0;

for i = 1:size(fileNames)
   currentSize = fileNames(i).bytes;
   if strfind(fileNames(i).name, '000.bin')
       currentSize = currentSize - headerSize;
   end
   totalSize = totalSize + currentSize;
end

totalSamples = totalSize/sampleSize;
totalTime = totalSamples/20000;
