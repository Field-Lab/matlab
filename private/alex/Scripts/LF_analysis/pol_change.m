function pol_change(ifON)

global fitParams
fitParams=calcParams(currentLFs,SpikeCount,handles);
if ifON % on
    [fitParams.amp, ind]=max(currentLFs(51:200,:)); 
    fitParams.pol=1;
    fitParams.currentPolaruty=1;
else % off
    [fitParams.amp, ind]=min(currentLFs(51:200,:)); 
    fitParams.pol=0;
    fitParams.currentPolaruty=-1;
end
fitParams.ind=ind+50;
fitParams.rightBorder=fitParams.ind+50;
fitParams.leftBorder=fitParams.ind-50;

ishandle(sbpl(i))

if ishandle(fitLines); delete(fitLines); end
for i=1:length(sbpl)
    realTrial=start_trial+i-1;
    subplot(sbpl(i));
    a=get(gca,'YLim');
    fitLines(2,i)=line([fitParams.ind(realTrial),fitParams.ind(realTrial)],a,'color','k','linewidth',2);
    fitLines(1,i)=line([fitParams.rightBorder(realTrial),fitParams.rightBorder(realTrial)],a,'color','g','linewidth',2);
    fitLines(3,i)=line([fitParams.leftBorder(realTrial),fitParams.leftBorder(realTrial)],a,'color','g','linewidth',2);
end