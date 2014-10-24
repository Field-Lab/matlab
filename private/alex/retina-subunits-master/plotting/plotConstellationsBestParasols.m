clear all
path = '~/Dropbox/Projects/retinaSubunits/code3/examples';
clf
set(gcf,'Position',[552 83 722 723]);


subp(2,1,1,0.01);
dataSet = '/2011-07-05-2/data002/';
file = '3203-off-parasol.mat';
load(fullfile(path,dataSet,file));
plotPairwiseCell(cellDat,cellDat.v,0)

subp(2,1,2,0.01);
dataSet = '/2010-03-05-2/data013/';
file = '211-on-parasol.mat';
load(fullfile(path,dataSet,file));
plotPairwiseCell(cellDat,cellDat.v,1)
