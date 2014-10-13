function splitHeaderData(dataFolderIn,dataFolderOut,visionPath)
% This function separates the data in a bin file from its header in order
% to facilitate bin files recombining.
% 
% SPLITHEADERDATA(dataFolderIn,dataFolderOut) will read the bin file in
% dataFolderIn and copy it to dataFolderOut. Two bin files will be created
% in DataFolderOut: DataFolderOut000.bin containing the header file only
% and DataFolderOut001.bin containing the data only.
%
% SPLITHEADERDATA(...,visionPath) will try to link to the vision jar file
% specified in VisionPath. If not specified the vision jar archive will be
% read from a default location.
% 
% Version: v4.04 - 05/08/2011

if nargin==2
    visionPath = '';    % Add here a default Vision path to not have to specify it every time
end

if ~exist('edu/ucsc/neurobiology/vision/io/ModifyRawDataFile','class')
    javaaddpath(visionPath);
end

% Linking to the relevant files
oldFile = edu.ucsc.neurobiology.vision.io.RawDataFile(dataFolderIn);
header = oldFile.getHeader();
if ~exist(dataFolderOut,'dir')
    mkdir(dataFolderOut);
end
newFile = edu.ucsc.neurobiology.vision.io.ModifyRawDataFile(dataFolderOut, header);
newFile.addFile();

% Copying the data to the new file by chunks of 10000 samples
nSamplesPerPulse = 10000;
nSamplesTotal = header.getNumberOfSamples();
startSample = 0;

while (startSample + nSamplesPerPulse) < nSamplesTotal
    data = oldFile.getData(startSample,nSamplesPerPulse);
    newFile.appendDataToLastFile(data);
    
    startSample = startSample + nSamplesPerPulse;
end

nSamplesLastPulse = nSamplesTotal - startSample;
if nSamplesLastPulse > 0
    data = oldFile.getData(startSample,nSamplesLastPulse);
    newFile.appendDataToLastFile(data);
end

% Closing the files
oldFile.close();
newFile.close();

end % splitHeaderData