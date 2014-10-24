function [train test celldat errFlag] = loadCellData(dat,celldat,plotting)

% dat is a data structure containing info for the full data set
% and for the cell we want to load
%
% loadType controls two methods of loading cells
% 1 - index into the desired cell type, e.g. the 1st on parasol, useful for
%     batch processing
% 2 - index into the raw cell id (e.g. from Vision)

errFlag = 0;

% copy stuff from cell data structure that we might need
celldat.dataSet = dat.dataSet;
celldat.dataPath = dat.dataPath;
celldat.analPath = dat.analPath;
celldat.figurePath = dat.figurePath;
celldat.mainDataPath = dat.mainDataPath;
celldat.mainAnalPath = dat.mainAnalPath;
celldat.mainFigurePath = dat.mainFigurePath;
celldat.computer = dat.computer;
celldat.locations = dat.locations;

if celldat.loadType == 1
  ids = getTypeIds(dat.cellTypes,celldat.cellType);
  cellInd = getCellInds(dat.rgcIds,ids(celldat.cellNum));
  celldat.rgcId = ids(celldat.cellNum);
end
if celldat.loadType == 2
  cellInd = find(dat.rgcIds == celldat.rgcId);
end
celldat.cellInd = cellInd;

if isempty(dat.stas.time_courses{celldat.cellInd})
  fprintf('(loadCellData) missing time course for cell %g\n',celldat.rgcId);
  train = [];
  test = [];
  errFlag = 1;
  return
end
celldat.cellType = getCellType(dat.cellTypes,dat.rgcIds(celldat.cellInd));
celldat.cellType = lower(celldat.cellType{1});
timeCourse = dat.stas.time_courses{celldat.cellInd};
polarity = dat.stas.polarities{celldat.cellInd};
if polarity == -1
  timeCourse = -timeCourse;
end
if size(timeCourse,2) > 1
  timeCourse = mean(timeCourse,2);
end
celldat.timeCourse = timeCourse;
celldat.polarity = polarity;

% get spike data
spikeRate = dat.spikeRate(celldat.cellInd,:);

if ~celldat.tempFilter
  spikeRate = circshift(spikeRate,[0 -1]);
end
spikeRate = spikeRate(:);

% get cone ids for the cell based on the threshold
if isfield(dat,'coneWeights')
  if sum(abs(dat.coneWeights(:,celldat.cellInd))) == 0
    fprintf('(loadCellData) missing cone weights for cell %g\n',celldat.rgcId);
    train = [];
    test = [];
    celldat = [];
    errFlag = 1;
    return
  end 
  celldat = getConeInds(dat,celldat);   
end

if size(celldat.locs_c,1)==0
  fprintf('(loadCellData) no cones for cell %g\n',celldat.rgcId);
  train = [];
  test = [];
  celldat = [];
  errFlag = 1;
  return;
end

if isfield(dat,'coneTypes');
  celldat.coneTypesFit = dat.coneTypes(celldat.keepInds);
end

celldat.staFit = celldat.sta(celldat.keepInds);
celldat.snr = max(abs(celldat.sta))/(1.4826*mad(celldat.sta,1));
celldat.rf = dat.rfs{celldat.cellInd};
celldat.xRange = ([floor(min(celldat.locs_c(:,1))-1.5) ...
  ceil(max(celldat.locs_c(:,1))+1.5)]);
celldat.yRange = ([floor(min(celldat.locs_c(:,2))-1.5) ...
  ceil(max(celldat.locs_c(:,2))+1.5)]);

celldat.xRange(1) = max(celldat.xRange(1),ceil(min(dat.locations(:,1))));
celldat.xRange(2) = min(celldat.xRange(2),floor(max(dat.locations(:,1))));
celldat.yRange(1) = max(celldat.yRange(1),ceil(min(dat.locations(:,2))));
celldat.yRange(2) = min(celldat.yRange(2),floor(max(dat.locations(:,2))));

if ~isempty(celldat.rf)
  celldat.rf = celldat.rf(celldat.yRange(1):celldat.yRange(2),...
    celldat.xRange(1):celldat.xRange(2),:);
end

% get training and testing inputs and response
if celldat.zscore
  coneInputs = zscore(dat.coneInputs(:,celldat.keepInds));
else
  coneInputs = dat.coneInputs(:,celldat.keepInds);
end

if celldat.tempFilter
  coneInputs = tempFilterCones(coneInputs,timeCourse/norm(timeCourse));
end

[train.X_ct test.X_ct] = partitionDat(coneInputs...
  ,celldat.trainFrac,celldat.subDivide,1);
[train.R_t test.R_t] = partitionDat(spikeRate,...
  celldat.trainFrac,celldat.subDivide,1);

train.X_ct = train.X_ct';
test.X_ct = test.X_ct';
train.R_t = train.R_t';
test.R_t = test.R_t';

if size(train.X_ct,2) == 0
  fprintf('(loadCellData) no valid cones for this cell\n');
elseif plotting
  figure
  set(gcf,'Position',[0 433 400 1000]);
  set(gcf,'Color','white');
  colormap(cjet);
  set(gca,'YDir','reverse');
  
  subplot(3,1,1);
  if isfield(dat,'vornoiMap')
    hold on
    image(dat.vornoiMap);
    set(gca,'YDir','reverse');
  end
  
  scatter(dat.locations(:,1),dat.locations(:,2),20,celldat.sta,'filled');
  axis equal
  xlim([min(dat.locations(:,1))*0.95 max(dat.locations(:,1))*1.05]);
  ylim([min(dat.locations(:,2))*0.95 max(dat.locations(:,2))*1.05]);
  set(gca,'CLim',[-max(abs(celldat.sta)) max(abs(celldat.sta))]);
  set(gca,'YDir','reverse')
  
  subplot(3,1,2);
  if isfield(dat,'vornoiMap')
    hold on
    image(dat.vornoiMap);
    set(gca,'YDir','reverse');
  end
  scatter(celldat.locs_c(:,1),celldat.locs_c(:,2),25,celldat.staFit,'filled');
  axis equal
  xlim([min(celldat.locs_c(:,1))*0.95 max(celldat.locs_c(:,1))*1.05]);
  ylim([min(celldat.locs_c(:,2))*0.95 max(celldat.locs_c(:,2))*1.05]);
  set(gca,'CLim',[-max(abs(celldat.staFit)) max(abs(celldat.staFit))]);
  if isfield(dat,'cellType')
    title(sprintf('cell# %g, %s',celldat.cellNum,celldat.cellType));
  end
  set(gca,'YDir','reverse')
  
  if exist('timeCourse','var');
    subplot(3,1,3);
    plot(timeCourse,'.-');
  end
  
end