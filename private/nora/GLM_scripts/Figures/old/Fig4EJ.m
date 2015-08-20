load('mat/1471rk1.mat')

plotraster_nolabel(FG.cross,FG.NSEM);
pause()
plotraster_nolabel(FG.NSEM.xvalperformance,FG.NSEM)
pause()
plotraster_nolabel(FG.WN.xvalperformance,FG.WN)
pause()

figure(1);
time=1/120*1000*(0:29);
plot(time,FG.WN.linearfilters.Stimulus.time_rk1,'k','LineWidth',1.5)
hold on
plot([0, 250],[0 0],'k')
xlabel('Time (ms)')

pause()

figure(1);
time=1/120*1000*(0:29);
plot(time,FG.NSEM.linearfilters.Stimulus.time_rk1,'k','LineWidth',1.5)
hold on
plot([0, 250],[0 0],'k')
xlabel('Time (ms)')

pause()
plotfilter(FG.WN,'PS');

pause() 
plotfilter(FG.NSEM,'PS');