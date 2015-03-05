function calcParams(currentLFs,SpikeCount,handles,trialsPerScreen,start_trial)

global fitParams


if get(handles.on,'Value') %on cell
    [amp,ind]=max(currentLFs(51:200,:));
    fitParams.currentPolarity=1;
    pol=1;
else
    [amp,ind]=min(currentLFs(51:200,:));
    fitParams.currentPolarity=-1;
    pol=0;
end
fitParams.ind=ind+50;
fitParams.rightBorder=fitParams.ind+50;
fitParams.leftBorder=fitParams.ind-50;
fitParams.pol=pol;
fitParams.amp=amp;
fitParams.SpikeCount=SpikeCount;

fitParams.limits=zeros(3,2);
for i=1:6
    fitParams.limits(i)=str2num(get(handles.getFitLimits(i),'String'));
end
fitParams.limits=fitParams.limits';

makePlots(currentLFs,trialsPerScreen,start_trial,[],1);