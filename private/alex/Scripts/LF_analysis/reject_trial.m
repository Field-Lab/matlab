function reject_trial(handles,fitParams,currentUnit,currentLFs)

global fitResults

a=str2num(get(handles.from,'String'));
b=str2num(get(handles.to,'String'));
if a==0 || b==0
    ginput(1);    
    find(fitParams.subplots(:,1)==gca)    
    trial=fitParams.subplots(fitParams.subplots(:,1)==gca,2);
    fitResults(:,trial,currentUnit)=[0 0 0];    
    subplotsToChange=find(fitParams.subplots(:,1)==gca);    
else
    trial=fitParams.subplots(a:b,2);
    fitResults(:,trial,currentUnit)=fitResults(:,fitParams.subplots(a:b,2),currentUnit)*0;
    subplotsToChange=a:b;
end

cnt=1;
for i=subplotsToChange
    subplot(fitParams.rows,fitParams.cols,i)    
    hold off
    plot(currentLFs(:,i));
    axis tight
    tmp=get(gca,'YLim');
    line([fitParams.ind(trial(cnt)),fitParams.ind(trial(cnt))],[tmp(1),tmp(2)],'color','k','linewidth',2);
    line([fitParams.rightBorder(trial(cnt)),fitParams.rightBorder(trial(cnt))],[tmp(1),tmp(2)],'color','g','linewidth',2);
    line([fitParams.leftBorder(trial(cnt)),fitParams.leftBorder(trial(cnt))],[tmp(1),tmp(2)],'color','g','linewidth',2);
    title(['trial ', int2str(trial(cnt)),' (',int2str(i),'), ',int2str(fitParams.SpikeCount(trial(cnt))),' sp']);

end