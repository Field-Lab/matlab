a=dir('S:\data\alexandra\MEA_data\20121216\Filters\*.mat')
i=1;
nd=1;
nds='8765432123333';
hc=0;
figure
set(gcf,'position',[4    34   810   939])
pl=1;
while i<=400
    if i>1
        if pl==1
            pl1=48;
        else
            pl1=pl-1;
        end
        subplot(8,6,pl1)             
        plot(mean(filters(:,1:60),2)/10e3,'color','b')
        line([0 500],[0,0],'color','k')
        title(['ND',p,' ',contr,' ',int2str(hc)])
        axis tight
        set(gca,'xtick',0,'xticklabel','')
    end
    if hc==12&&strcmp(contr,'Low')
        hc=0;
    end
    load(['S:\data\alexandra\MEA_data\20121216\Filters\',a(i).name])
    if mod(i,2)
        contr='High';
        hc=hc+1;
    else
        contr='Low';
    end
    
    subplot(8,6,pl)
    plot(mean(filters(:,1:60),2)/10e3,'color','r','linewidth',2)
    line([0 500],[0,0],'color','k')
    axis tight
    p=nds(nd);
    title(['ND',p,' ',contr,' ',int2str(hc)])
    set(gca,'xtick',0,'xticklabel','')
    
    if mod(i,24)==0
        nd=nd+1;
    end    
    
    i=i+1;
    while i>=length(a)+1
        pause(15)
        a=dir('S:\data\alexandra\MEA_data\20121216\Filters\*.mat');
        length(a)
    end
    pl=pl+1;
    if pl==49
        pl=1;
    end
end


% 
% for i=1:60
%     subplot(6,10,i)
%     plot(filters(:,i))
% end
