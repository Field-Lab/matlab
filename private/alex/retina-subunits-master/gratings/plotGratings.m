function plotGratings(g,grat,dat,mode,gFit,newFig)

%-------------------------------------------
%
% g = plotGratings(g,grat,dat,mode)
%
% plot the output of a grating analysis for a single cell
%
% g is a structure of grating response parameters
% for a single cell (see 'getGratingResp.m')
%
% grat needs to contain the fields:
% grat.stimulus -- structure with stimulus info
%
% dat needs to contain the fields:
% dat.locs_c -- cone xy locations for this cell, cones x 2
% dat.xRange -- range of x values including the cell's cones
% dat.yRange -- range of y values including the cell's cones
%
% mode has three options:
% 'cycle' -- plot cycle avearged responses
% 'raster' -- plot a raster for each cycle
% 'freq' -- plot the f1 and f2 versus spatial frequency
%
% freeman, 5-5-2012
%-------------------------------------------

% set some constants
REFRESH = 60;
STIMDUR = 8;

% get some parameters
nPh = size(g.cycleAve,3);
nPer = size(g.cycleAve,2);
nReps = grat.stimulus.repetitions;
cycleDur = grat.stimulus.params.TEMPORAL_PERIOD;
nCycles = (STIMDUR*REFRESH)/cycleDur;

% are we going to plot the cones?
if exist('dat','var') && ~isempty(dat)
  plotCones = 1;
else
  plotCones = 0;
end

if ~exist('newFig','var') || isempty(newFig);
  newFig = 1;
end

% figure out what kind of plot we're making
switch mode
  case 'cycle'
    if newFig
      clf;
      set(gcf,'Color','white');
      set(gcf,'Position',[1 66 1440 740]);
      colormap(gray);
    end
    mx = round2(max(g.cycleAve(:)),0.1,'ceil');
    disp(mx)
    for iph=1:nPh
      for iper=1:nPer
        subp(nPh*2,nPer,nPer*(iph*2-1)+iper,0.01);
        % was shifting by -19 for 2010-03-05-2/data016
        % was shifting by -9 for 2011-10-25-5/data001
        h = plot(circshift(squeeze(g.cycleAve(:,iper,iph)),[-9 0]),'k','LineWidth',3);
        if exist('gFit','var');
          hold on
          hFit = plot(circshift(squeeze(gFit.cycleAve(:,iper,iph)),[-9 0]),'r','LineWidth',3);
        end
        ylim([0 mx*0.6]); % was using 0.6
        xlim([1 15]);
        set(gca,'TickDir','out');
        set(gca,'XTick',[1 30]);
        set(gca,'YTick',[0 mx]);
        if ~newFig
          set(h,'Color','k');
        end
        box off;
        if ~(iper==1 && iph==nPh)
          set(gca,'XTickLabel',[]);
          set(gca,'YTickLabel',[]);
        else
          set(gca,'XTickLabel',[0 0.25]);
          set(gca,'YTickLabel',[0 mx]);
        end
        if plotCones && newFig
          subp(nPh*2,nPer,2*nPer*(iph-1)+iper,0.01);
          scatter(dat.locs_c(:,1),dat.locs_c(:,2),100,-1*g.coneMaxFrame(:,iper,iph),'filled');
          set(gca,'CLim',[-1 1]);
          set(gca,'YDir','reverse');
          axis equal
          xlim(dat.xRange); ylim(dat.yRange);
          axis off
        end
      end
    end
    
  case 'raster'
    clf;
    set(gcf,'Color','white');
    set(gcf,'Position',[1 66 1440 740]);
    colormap(redgrayblue);
    for iph=1:nPh
      for iper=1:nPer
        subp(nPh*2,nPer,nPer*(iph*2-1)+iper,0.01);
        raster = reshape(g.timeCourses(:,:,iper,iph),[cycleDur nCycles*nReps]);
        raster = clip(raster,0,2);
        flipped = 2-raster';
        scaled = flipped/max(flipped(:));
        rasterPlot(:,:,1) = scaled;
        rasterPlot(:,:,2) = scaled;
        rasterPlot(:,:,3) = scaled;
        image(rasterPlot);
        xlim([1 30]);
        set(gca,'YDir','reverse');
        box off
        if ~(iper==1 && iph==nPh)
          set(gca,'TickDir','out');
          set(gca,'XTick',[]);
          set(gca,'YTick',[]);
        else
          set(gca,'TickDir','out');
          set(gca,'XTick',[1 30]);
          set(gca,'XTickLabel',[0 0.25]);
          set(gca,'YTick',[]);
          set(gca,'YTickLabel',[]);
        end
        if plotCones
          subp(nPh*2,nPer,2*nPer*(iph-1)+iper,0.01);
          scatter(dat.locs_c(:,1),dat.locs_c(:,2),20,g.coneMaxFrame(:,iper,iph),'filled');
          set(gca,'CLim',[-1 1]);
          set(gca,'YDir','reverse');
          axis equal
          xlim(dat.xRange); ylim(dat.yRange);
          axis off
        end
      end
    end
    
  case 'freq'
    clf;
    set(gcf,'Color','white');
    set(gcf,'Position',[900 465 475 350]);
    hold on
    plot((1000/5.8)./fliplr(grat.stimulus.params.SPATIAL_PERIOD),mean(squeeze(mean(g.f1,1)),2),'k.-');
    plot((1000/5.8)./fliplr(grat.stimulus.params.SPATIAL_PERIOD),mean(squeeze(mean(g.f2,1)),2),'r.-');
    set(gca,'XScale','log');
    ylim([0 round2(max([g.f1(:);g.f2(:)]),0.1,'ceil')]);
end


