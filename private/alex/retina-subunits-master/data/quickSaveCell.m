% script for saving out a data file that contains
% everything we need to analyze a single cell
clear all
dat = loadData('laptop','2010-03-05-2/data013');

%%
path = '~/Dropbox/Projects/retinaSubunits/code3/examples';
mkdir(path,dat.dataSet);

%%
for i=[677];
dat.cellType = 'off midget';
dat.cellNum = i; %4621 and 2822 for 12-13-2011
dat.trainFrac = 0.8;
dat.subDivide = 10;
dat.thresh = 2;
dat.percent = 0.2;
dat.threshType = 'peak';
dat.zscore = 1;
dat.tempFilter = 1;
[train test dat] = loadCellData(dat,2,0);
v = loadVmat(dat,'vMats');

cellDat.dataPath = dat.dataPath;
cellDat.analPath = dat.analPath;
cellDat.dataSet = dat.dataSet;
cellDat.cellType = dat.cellType{1};
cellDat.rgcId = dat.rgcId;
cellDat.trainFrac = dat.trainFrac;
cellDat.subDivide = dat.subDivide;
cellDat.percent = dat.percent;
cellDat.threshType = dat.threshType;
cellDat.zscore = dat.zscore;
cellDat.tempFilter = dat.tempFilter;
cellDat.locs_c = dat.locs_c;
cellDat.coneIds = dat.conesFit;
cellDat.center = dat.center;
cellDat.sta = dat.staFit;
cellDat.v = v;

rf = dat.rfs{dat.cellInd};
xRange = round([min(cellDat.locs_c(:,1))*0.975 max(cellDat.locs_c(:,1))*1.025]);
yRange = round([min(cellDat.locs_c(:,2))*0.975 max(cellDat.locs_c(:,2))*1.025]);
if ~isempty(rf)
    rf = rf(yRange(1):yRange(2),xRange(1):xRange(2),:);
end
cellDat.rf = rf;

saveName = sprintf('%g-%s.mat',cellDat.rgcId,strrep(lower(cellDat.cellType),' ','-'));
save(fullfile(path,dat.dataSet,saveName),'train','test','cellDat');
end