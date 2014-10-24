function batchPlot(computer,dataSet,cellType,modelType,percent)

dat = loadData(computer,dataSet);

%% plot the cones for this piece
clf
switch dat.computer
  case 'laptop'
    set(gcf,'Position',[1 1 620 800]);
  case 'lcv'
    set(gcf,'Position',[1 77 770 1020]);
end
set(gcf,'Position',[1 1 620 800]);
h1 = axes('Position',[0.01 0.01 0.98 0.98]);
xPos = dat.locations(:,1);
yPos = dat.locations(:,2);
xRangeFull = [min(xPos),max(xPos)];
xRangeFull = xRangeFull + [-0.02*diff(xRangeFull) 0.02*diff(xRangeFull)];
yRangeFull = [min(yPos),max(yPos)];
yRangeFull = yRangeFull + [-0.02*diff(yRangeFull) 0.02*diff(yRangeFull)];
rect1 = rectangle('Position',[xRangeFull(1) yRangeFull(1) diff(xRangeFull) diff(yRangeFull)]);
set(rect1,'FaceColor',[0.5 0.5 0.5]);
set(rect1,'EdgeColor','none');
set(gca,'YDir','reverse');
for i=1:length(xPos)
  circle(xPos(i),yPos(i),1,'k','none',0.5);
end
axis equal;
axis off;
xlim(xRangeFull);
ylim(yRangeFull);
%% plot the cells
numCells = getCellTypeNum(dat.cellTypes,cellType);
for i=1:numCells
  %try
    fprintf('(batchPlot) analyzing %s %g/%g\r',cellType,i,numCells);
    
    celldat.cellType = cellType;
    celldat.cellNum = i;
    celldat.loadType = 1;
    celldat.trainFrac = 0.8;
    celldat.subDivide = 10;
    celldat.thresh = 15;
    celldat.percent = percent;
    celldat.threshType = 'peak';
    celldat.zscore = 1;
    celldat.tempFilter = 1;
    [train test celldat] = loadCellData(dat,celldat,0);
    if ~isempty(celldat)
      res = loadFit(celldat,modelType,[]);
      if ~isempty(res)
        if res.out_SUB.r2 >= 0.05
          plotSubunitsLinesOvals(res.fit_SUB,celldat,1.4,2.2,0.3);
        end
      end
    end
    drawnow;
  %catch
  %  fprintf('(batchPlot) skipping %s %g/%g\r',cellType,i,numCells);
  %end
end
%%
set(gcf,'InvertHardcopy', 'off');
set(gcf,'PaperPositionMode', 'auto');
saveName = strcat(dat.dataSet(1:12),'-',dat.dataSet(14:20),'-',cellType);
saveName = strcat(saveName,'-',num2str(percent));
print('-dpdf',fullfile(dat.mainPath,'populations',strcat(saveName,'.pdf')));
