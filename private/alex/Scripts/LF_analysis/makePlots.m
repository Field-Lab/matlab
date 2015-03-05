function makePlots(currentLFs,trialsPerScreen,start_trial,resFit,flag)

global fitParams

end_trial=min(start_trial+trialsPerScreen-1,size(currentLFs,2));

[rows,cols]=opt_subplots(trialsPerScreen);
fitParams.rows=rows;
fitParams.cols=cols;
cnt=1;
for trial=start_trial:end_trial
    
    fitParams.subplots(cnt,1)=subplot(rows,cols,cnt);
    fitParams.subplots(cnt,2)=trial;
    if flag==1 % replot completely
        hold off
        plot(currentLFs(:,trial));
        axis tight
        tmp=get(gca,'YLim');
        
        line([fitParams.ind(trial),fitParams.ind(trial)],[tmp(1),tmp(2)],'color','k','linewidth',2);
        line([fitParams.rightBorder(trial),fitParams.rightBorder(trial)],[tmp(1),tmp(2)],'color','g','linewidth',2);
        line([fitParams.leftBorder(trial),fitParams.leftBorder(trial)],[tmp(1),tmp(2)],'color','g','linewidth',2);
        title(['trial ', int2str(trial),' (',int2str(cnt),'), ',int2str(fitParams.SpikeCount(trial)),' sp']);
    else % plot fit results
        hold on 

        x=(fitParams.leftBorder(trial):fitParams.rightBorder(trial))';
        y=resFit(1,cnt)*exp(-((x-resFit(2,cnt))/resFit(3,cnt)).^2);
        plot(x,y,'r','LineWidth',2)
        a=int2str(round(resFit(:,cnt))');
        a(regexp(a,'  '))='';
        a(regexp(a,'  '))='';
        text(200,0.8*fitParams.amp(trial),a,'FontSize',10);
 
    end
    cnt=cnt+1;
end


