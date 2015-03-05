function change_id(currentLFs,flag)
global fitParams


if flag
    [tmp,~]=ginput(1);
    tmp=round(tmp);
    trial=fitParams.subplots(fitParams.subplots(:,1)==gca,2);
    
    fitParams.ind(trial)=tmp;
    fitParams.rightBorder(trial)=min(tmp+50,430);
    fitParams.leftBorder(trial)=max(30,tmp-50);
else    
    [tmp,~]=ginput(2);
    tmp=round(tmp);
    trial=fitParams.subplots(fitParams.subplots(:,1)==gca,2);
    
    fitParams.leftBorder(trial)=max(30,min(tmp));
    fitParams.rightBorder(trial)=min(430,max(tmp));
   
end

subplot(fitParams.rows,fitParams.cols,find(fitParams.subplots(:,1)==gca))
hold off
plot(currentLFs(:,trial));
axis tight
tmp=get(gca,'YLim');
line([fitParams.ind(trial),fitParams.ind(trial)],[tmp(1),tmp(2)],'color','k','linewidth',2);
line([fitParams.rightBorder(trial),fitParams.rightBorder(trial)],[tmp(1),tmp(2)],'color','g','linewidth',2);
line([fitParams.leftBorder(trial),fitParams.leftBorder(trial)],[tmp(1),tmp(2)],'color','g','linewidth',2);
title(['trial ', int2str(trial),' (',int2str(find(fitParams.subplots(:,1)==gca)),'), ',int2str(fitParams.SpikeCount(trial)),' sp']);
