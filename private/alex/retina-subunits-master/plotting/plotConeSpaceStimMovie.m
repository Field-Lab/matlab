clf
set(gcf,'Position',[0 0 1400 510]);
set(gcf,'Color','w');
%%
vidObj = VideoWriter('stimmovie.avi');
vidObj.FrameRate = 2;
vidObj.Quality = 100;
%%
open(vidObj);
for i=1:20
  
  rndframe = (rand(16,15)-0.5)*0.8 + 0.5;
  rndframe = ones(19,15);
  
  %% plot the rf in pixel space
  
  coneSz = 1.4;
  subSz = 2;
  coneEdge = 1.3;
  figWidth = 0.3;
  figHeight = 0.9;
  
  % plot the RF in pixel space
  if isfield(celldat,'rf') && ~isempty(celldat.rf)
    h = axes('Position',[0.35 0.05 figWidth figHeight]);
    image(norm_image(mean(celldat.rf,3)));
    xlim(celldat.xRange - celldat.xRange(1));
    ylim(celldat.yRange - celldat.yRange(1));
    axis image
    set(h,'YDir','reverse');
    axis off
  end
  %%
  % plot the cones for this cell
  hold on
  xPos = celldat.locs_c(:,1) - min(celldat.xRange) + 1;
  yPos = celldat.locs_c(:,2) - min(celldat.yRange) + 1;
  for i=1:length(xPos)
    circle(xPos(i),yPos(i),coneSz,'none','k',coneEdge);
  end
  
  
  %%
  h = axes('Position',[0.025 0.05 figWidth figHeight]);
  image(repmat(rndframe,[1 1 3]));
  xlim(celldat.xRange - celldat.xRange(1));
  ylim(celldat.yRange - celldat.yRange(1));
  axis image off
  set(h,'YDir','reverse');
  
  %%
  h2 = axes('Position',[0.675 0.05 figWidth figHeight]);
  
  axis equal
  xlim(celldat.xRange);
  ylim(celldat.yRange);
  set(gca,'YDir','reverse');
  axis off
  
  rect1 = rectangle('Position',[celldat.xRange(1) celldat.yRange(1) sum(celldat.xRange) sum(celldat.yRange)]);
  set(rect1,'FaceColor',[0.5 0.5 0.5]);
  set(rect1,'EdgeColor','none');
  
  
  % plot the cones
  hold on
  xPos = celldat.locs_c(:,1);
  yPos = celldat.locs_c(:,2);
  
  coneinds(:,1) = round(celldat.locs_c(:,1)-celldat.xRange(1));
  coneinds(:,2) = round(celldat.locs_c(:,2)-celldat.yRange(1));
  
  for i=1:size(coneinds,1)
    coneresponses(i) = rndframe(coneinds(i,1),coneinds(i,2));
  end
  
  for i=1:size(celldat.locs_c,1)
    clr = repmat(coneresponses(i),1,3);;
    circle(xPos(i),yPos(i),coneSz,clr,'k',coneEdge);
  end
  drawnow;
  %%
  currFrame = getframe(gcf);
  %%
  writeVideo(vidObj,currFrame);
end
close(vidObj);