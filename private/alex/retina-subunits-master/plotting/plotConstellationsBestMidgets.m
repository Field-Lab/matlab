clear all
path = '~/Dropbox/Projects/retinaSubunits/code3/examples';
clf
set(gcf,'Position',[552 83 722 723]);

subp(4,4,1,0.01);
dataSet = '/2010-03-05-2/data013/';
file = '677-on-midget.mat';
load(fullfile(path,dataSet,file));
plotPairwiseCell(cellDat,cellDat.v,1)

subp(4,4,2,0.01);
dataSet = '/2010-03-05-2/data013/';
file = '1247-on-midget.mat';
load(fullfile(path,dataSet,file));
plotPairwiseCell(cellDat,cellDat.v,1)

subp(4,4,3,0.01);
dataSet = '/2010-03-05-2/data013/';
file = '934-on-midget.mat';
load(fullfile(path,dataSet,file));
plotPairwiseCell(cellDat,cellDat.v,1)

subp(4,4,4,0.01);
dataSet = '/2010-03-05-2/data013/';
file = '1085-on-midget.mat';
load(fullfile(path,dataSet,file));
plotPairwiseCell(cellDat,cellDat.v,1)

subp(4,4,5,0.01);
dataSet = '/2010-03-05-2/data013/';
file = '1156-off-midget.mat';
load(fullfile(path,dataSet,file));
plotPairwiseCell(cellDat,cellDat.v,0)

subp(4,4,6,0.01);
dataSet = '/2010-03-05-2/data013/';
file = '1352-on-midget.mat';
load(fullfile(path,dataSet,file));
plotPairwiseCell(cellDat,cellDat.v,1)

subp(4,4,7,0.01);
dataSet = '/2008-08-27-5/data003/';
file = '1112-off-midget.mat';
load(fullfile(path,dataSet,file));
plotPairwiseCell(cellDat,cellDat.v,0)

subp(4,4,8,0.01);
dataSet = '/2010-03-05-2/data013/';
file = '1516-on-midget.mat';
load(fullfile(path,dataSet,file));
plotPairwiseCell(cellDat,cellDat.v,1)

subp(4,4,9,0.01);
dataSet = '/2010-03-05-2/data013/';
file = '1636-on-midget.mat';
load(fullfile(path,dataSet,file));
plotPairwiseCell(cellDat,cellDat.v,1)

subp(4,4,10,0.01);
dataSet = '/2010-03-05-2/data013/';
file = '4936-off-midget.mat';
load(fullfile(path,dataSet,file));
plotPairwiseCell(cellDat,cellDat.v,0)

subp(4,4,11,0.01);
dataSet = '/2008-08-27-5/data003/';
file = '1411-off-midget.mat';
load(fullfile(path,dataSet,file));
plotPairwiseCell(cellDat,cellDat.v,0)

subp(4,4,12,0.01);
dataSet = '/2011-07-05-2/data002/';
file = '2837-on-midget.mat';
load(fullfile(path,dataSet,file));
plotPairwiseCell(cellDat,cellDat.v,1)

subp(4,4,13,0.01);
dataSet = '/2008-08-27-5/data003/';
file = '528-off-midget.mat';
load(fullfile(path,dataSet,file));
plotPairwiseCell(cellDat,cellDat.v,0);

subp(4,4,14,0.01);
dataSet = '/2008-08-27-5/data003/';
file = '542-off-midget.mat';
load(fullfile(path,dataSet,file));
plotPairwiseCell(cellDat,cellDat.v,0)

subp(4,4,15,0.01);
dataSet = '/2008-08-27-5/data003/';
file = '691-off-midget.mat';
load(fullfile(path,dataSet,file));
plotPairwiseCell(cellDat,cellDat.v,0)

subp(4,4,16,0.01);
dataSet = '/2008-08-27-5/data003/';
file = '1052-off-midget.mat';
load(fullfile(path,dataSet,file));
plotPairwiseCell(cellDat,cellDat.v,0)

