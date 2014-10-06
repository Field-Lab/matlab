%javaaddpath 'J:\2010\Vision.jar'; %define path to Vision jar file

full_path='J:\2010\data005\data005000.bin'; %define path to raw data file
rawFile=edu.ucsc.neurobiology.vision.io.RawDataFile(full_path); 
RawData=rawFile.getData(0,20000)'; %the output is 513x20000 array. First index is channel number, and there is 40000 samples for each channel. The first sample is sample number 100000, as specified in the first argument.

channel=[3:6];
RawDataChannel=RawData(channel+1,:)';