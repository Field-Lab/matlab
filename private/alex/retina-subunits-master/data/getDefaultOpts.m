function celldat = getDefaultOpts(celldat)

% get a set of default options for loading a cell structure
% this stuff will rarely change

celldat.trainFrac = 0.8; % default by Jeremy
celldat.subDivide = 10; % default by Jeremy

% celldat.trainFrac = 0.7;
% celldat.subDivide = 20;

celldat.thresh = 15;
celldat.threshType = 'peak';
celldat.zscore = 1;
celldat.tempFilter = 1;