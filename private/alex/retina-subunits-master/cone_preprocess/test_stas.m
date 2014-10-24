clf
for i=1:getCellTypeNum(dat.cellTypes,'off midget');
  dat.cellType = 'off midget';
  dat.cellNum = i;
  dat.rgcId = 3; 
  dat.trainFrac = 0.8;
  dat.subDivide = 10;
  dat.thresh = 15;
  dat.percent = 0.33;
  dat.threshType = 'peak';
  dat.zscore = 1;
  dat.tempFilter = 1;
  [train test dat] = loadCellData(dat,1,0);

  mappedCell = datarun.cell_list_map(dat.cellInd);
  mappedCell
  if ~isempty(mappedCell{1})
    
    shift = -1;
    cell = find(datarun.cell_ids==mappedCell{1});
    sta = simpleSta(datarun.cone_inputs,circshift(datarun.spike_rate(cell,:),[0 shift])');
     scatter(dat.locations(:,1),dat.locations(:,2),20,sta,'filled');
     set(gca,'CLim',[-max(abs(sta(:))),max(abs(sta(:)))]);
     colormap(cjet);
     set(gca,'YDir','reverse');
     drawnow;
     
    weights = abs(sta)/sum(abs(sta));
    weightsTmp = weights;
    weightsTmp(weights<0.5*max(weights)) = 0;
    weightsTmp = weightsTmp/sum(weightsTmp);
    center = weightsTmp'*dat.locations;
    hold on
    plot(dat.center(1),dat.center(2),'k.','MarkerSize',20);
    plot(center(1),center(2),'r.','MarkerSize',20);
    plot([dat.center(1),center(1)],[dat.center(2),center(2)],'k-');
    drawnow;
    xlim([0 300]);
    ylim([0 300]);
  end
end

%%
set(gcf,'Position',[0 0 400 400]);
title('cones');
axis equal
xlim([0 300]);
ylim([0 300]);