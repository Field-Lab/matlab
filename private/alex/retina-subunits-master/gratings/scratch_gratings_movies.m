per = 9;
ph = 2;
clrs = gray(256);
set(gcf,'Position',[885 180 347 599]);
set(gcf,'Color','w');
%%
vidObj = VideoWriter('cone-low-freq.avi');
vidObj.FrameRate = 19; 
%%
open(vidObj);


ncycles = 1;
nframes = 30;
refPhase = 0;
x = [0:nframes-1];
sin_ts = sin(x*2*pi*(ncycles/nframes)-refPhase);
sin_ts = sin_ts(:);
sin_ts = circshift(sin_ts,[-16 0]);
y = circshift(g.cycleAve(:,per,ph),[-19 0]);
for i=1:30
  
  axes('Position',[0.1 0.55 0.8 0.4]);
  cla
  for j=1:celldat.nConesFit
    clrVal = g.coneMaxFrame(j,per,ph) * sin_ts(i);
    clr = clrs(ceil((clrVal + 1)/2 * 255),:);
    circle(celldat.locs_c(j,1),celldat.locs_c(j,2),1.5,clr,'k',1);
  end
  axis off
  axis equal
  set(gca,'YDir','reverse');
  xlim(celldat.xRange);
  ylim(celldat.yRange);
  
  %%
  subplot(2,1,2);
  cla
  plot(y(1:i),'Color',[0.2 0.2 0.2],'LineWidth',5);
  xlim([0 31]);
  ylim([0 0.55]);
  box off
  set(gca,'TickDir','out');
  set(gca,'FontSize',20);
  xlabel('Time (ms)');
  ylabel('Firing rate (ips)');
  set(gca,'YTick',[0 0.5]);
  
  currFrame = getframe(gcf);
  %%
  writeVideo(vidObj,currFrame);
  
end
close(vidObj);
