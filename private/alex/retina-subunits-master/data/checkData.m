function out = checkData(computer,dataSet)

switch computer
  case 'lcv'
    mainDataPath = '/share/users-lcv/freeman/retinalData/data/';
    mainAnalPath = '/share//users-lcv/freeman/retinalData/anal/';
    mainFigurePath = '/share/users-lcv/freeman/retinalData/figure/';
  case 'laptop'
    mainDataPath = '~/Dropbox/Projects/retinaSubunits/data/';
    mainAnalPath = '~/Dropbox/Projects/retinaSubunits/anal/';
    mainFigurePath = '~/Dropbox/Projects/retinaSubunits/figure/';
  case 'salk'
    mainDataPath = '/snle/acquisition/freeman/';
    mainAnalPath = '/Users/vision/Desktop/freeman/anal/';
    mainFigurePath = '/Users/vision/Desktop/freeman/figure/';
end

dataPath = strcat(mainDataPath,dataSet);
if ~isfile(fullfile(dataPath,'conepreprocess.mat'));
  error(sprintf('(checkData) data set %s missing',dataSet));
end

