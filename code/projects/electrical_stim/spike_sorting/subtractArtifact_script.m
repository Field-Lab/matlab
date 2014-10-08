labviewOutputFileFolder = '/Volumes/Data/2014-08-13-0/';
dataFolderOut= '/Volumes/Analysis/2014-08-13-0/data003-modified';
% dataFolderOut= '/Volumes/Analysis/2014-08-13-0/data012'; %Testing to see if vision input must be 7 characters long
dataFolderIn= '/Volumes/Data/2014-08-13-0/data003';
artDataFolderIn= '/Volumes/Data/2014-08-13-0/data010';
visionWritePath = '/Users/grosberg/code/write data GG/WriteDataFile.jar';

subtractArtRawData(dataFolderIn,artDataFolderIn, dataFolderOut,labviewOutputFileFolder, visionWritePath);

