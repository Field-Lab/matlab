unitsPath=dir(['S:\user\alexandra\MEA_data\',date,'\easy_formatted_units\*FFFlash*'])

nds='8765432123456';
for unit=1:length(unitsPath)
    load(['S:\user\alexandra\MEA_data\',date,'\easy_formatted_units\',unitsPath(unit).name])
    white_flash=zeros(55,4501);
    black_flash=zeros(55,4501);
    for trial=1:size(spike_info.flip_times,1)
        spikes=cell2mat(spike_info.spike_times(trial));
        flips=cell2mat(spike_info.flip_times(trial));
        flips=flips(:,1);
        conv=convolved(spikes,40,flips(end));
        conv=conv(121:end-120);
        for st=1:4:20
            white_flash(trial,1:4501)=white_flash(trial,1:4501)+conv(flips(st)-500:flips(st)+4000);
            black_flash(trial,1:4501)=black_flash(trial,:)+conv(flips(st+2)-500:flips(st+2)+4000);
        end
    end
    white_flash=white_flash/5;
    black_flash=black_flash/5;
    close(figure(1))
    figure(1)
    set(gcf,'position',[1  31 1680 946]);
    cnt=1;
    for nd=1:5:50
        subplot(5,4,cnt)
        plot(white_flash(nd:nd+4,:)','LineWidth',1.5);
        line([500,500],[0,max(max(white_flash(nd:nd+4,:)))],'Color','k','LineWidth',3);
        line([2500,2500],[0,max(max(white_flash(nd:nd+4,:)))],'Color','k','LineWidth',3);
        axis tight
        set(gca,'XTick',0)
        title(['ND',nds((nd+4)/5),' Flash RGB 50'],'FontSize',14,'FontWeight','bold','color',[0.6 0.6 0.6])
        if nd==1
            legend(('12345')')
        end
        subplot(5,4,cnt+1)
        plot(black_flash(nd:nd+4,:)','LineWidth',1.5);
        line([500,500],[0,max(max(black_flash(nd:nd+4,:)))],'Color','k','LineWidth',3);
        line([2500,2500],[0,max(max(black_flash(nd:nd+4,:)))],'Color','k','LineWidth',3);
        axis tight
        set(gca,'XTick',0)
        title(['ND',nds((nd+4)/5),' Flash RGB 10'],'FontSize',14,'FontWeight','bold')
        cnt=cnt+2;
    end
    subplot('Position',[0.5 0.95 0.00001 0.00001])
    set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
    title([unitsPath(unit).name(1:end-25),'. Spike Rate in response to FF flash 2s, BG30. Each line - mean of 5 flashes'],'FontSize',15,'FontWeight','bold','Interpreter','None')
    saveas(gcf,['F:\',date,'\FFFlash_plots\',unitsPath(unit).name(1:end-25),'_BG30.bmp']);

end




%% same, with equalized axis
unitsPath=dir(['S:\user\alexandra\MEA_data\',date,'\easy_formatted_units\*FFFlash*'])

nds='6543212';
for unit=1:length(unitsPath)
    load(['S:\user\alexandra\MEA_data\',date,'\easy_formatted_units\',unitsPath(unit).name])
    white_flash=zeros(50,4501);
    black_flash=zeros(50,4501);
    for trial=1:size(spike_info.flip_times,1)
        spikes=cell2mat(spike_info.spike_times(trial));
        flips=cell2mat(spike_info.flip_times(trial));
        flips=flips(:,1);
        conv=convolved(spikes,40,flips(end));
        conv=conv(121:end-120);
        for st=1:4:20
            white_flash(trial,1:4501)=white_flash(trial,1:4501)+conv(flips(st)-500:flips(st)+4000);
            black_flash(trial,1:4501)=black_flash(trial,:)+conv(flips(st+2)-500:flips(st+2)+4000);
        end
    end
    white_flash=white_flash/5;
    black_flash=black_flash/5;
    close(figure(1))
    figure(1)
    set(gcf,'position',[1  31 1680 946]);
    cnt=1;
    for nd=1:5:30
        subplot(3,4,cnt)
        plot(white_flash(nd:nd+4,:)','LineWidth',1.5);
        line([500,500],[0,max(max(white_flash))],'Color','k','LineWidth',3);
        line([2500,2500],[0,max(max(white_flash))],'Color','k','LineWidth',3);
        axis tight
        set(gca,'XTick',0)
        title(['ND',nds((nd+4)/5),' Flash RGB 50'],'FontSize',14,'FontWeight','bold','color',[0.6 0.6 0.6])
        if nd==1
            legend(('12345')')
        end
        subplot(3,4,cnt+1)
        plot(black_flash(nd:nd+4,:)','LineWidth',1.5);
        line([500,500],[0,max(max(black_flash))],'Color','k','LineWidth',3);
        line([2500,2500],[0,max(max(black_flash))],'Color','k','LineWidth',3);
        axis tight
        set(gca,'XTick',0)
        title(['ND',nds((nd+4)/5),' Flash RGB 10'],'FontSize',14,'FontWeight','bold')
        cnt=cnt+2;
    end
    subplot('Position',[0.5 0.95 0.00001 0.00001])
    set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
    title([unitsPath(unit).name(1:end-25),'. Spike Rate in response to FF flash 2s, BG30. Each line - mean of 5 flashes'],'FontSize',15,'FontWeight','bold','Interpreter','None')
    saveas(gcf,['F:\',date,'\FFFlash_plots_equalAxis\',unitsPath(unit).name(1:end-25),'_BG30.emf']);

end



%% old type of data - only long full field flashes of 2 contrast
date='20120405_1'
hekaPath=dir(['S:\user\alexandra\MEA_data\',date,'\HEKA\*phys'])
for i=1:2400
    a=hekaPath(i).name;
    b=regexp(a,'contrast_');
    d(i)=str2num(a(b+9));    
end
unitsPath=dir(['S:\user\alexandra\MEA_data\',date,'\easy_formatted_units\*FFFlash*'])
load(['S:\user\alexandra\MEA_data\',date,'\nd_by_files.mat'])
nd_change=[0 find(diff(c)~=0)];
nds='98765432123456787654321';
for unit=1:length(unitsPath)
    load(['S:\user\alexandra\MEA_data\',date,'\easy_formatted_units\',unitsPath(unit).name])
    flash=zeros(2400,4500);
    for trial=1:size(spike_info.flip_times,1)        
        spikes=cell2mat(spike_info.spike_times(trial));
        flips=cell2mat(spike_info.flip_times(trial));
        tmp=find(flips(:,2)~=-1)-2;
        flips=flips(:,1);
        conv=convolved(spikes,40,flips(end));
        flips=flips(tmp,1);
        conv=conv(121:end-120);
        flash(trial,1:4500)=conv(flips-499:flips+4000);
    end
    
   
    close(figure(1))
    figure(1)
    set(gcf,'position',[1  31 1680 946]);
    cnt=1;
    t=0;
    for nd=1:20        
        subplot(5,4,cnt)
        k=nd_change(nd)+1:2:nd_change(nd+1);
        toPlot=flash(k,:);
        clear toPlot1
        cnt1=1;
        for i=1:5:length(k)
            goto=min(i+4,size(toPlot,1));
            toPlot1(cnt1,1:4500)=mean(toPlot(i:goto,:));
            cnt1=cnt1+1;
        end
        plot(toPlot1','LineWidth',1.5);
        t=max(t,max(toPlot1(:)));
        set(gca,'XTick',0)
        title(['ND',nds(nd)],'FontSize',14,'FontWeight','bold')
        cnt=cnt+1;
    end
    for nd=1:20
        subplot(5,4,nd)
        line([500,500],[0,t],'Color','k','LineWidth',3);
        line([2500,2500],[0,t],'Color','k','LineWidth',3);
        axis tight
    end
    subplot('Position',[0.5 0.95 0.00001 0.00001])
    set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
    title([unitsPath(unit).name(1:end-25),'. Flash 40, BG 25. Each line - mean of 5 flashes'],'FontSize',15,'FontWeight','bold','Interpreter','None')
    saveas(gcf,['F:\',date,'\FFFlash_plots\',unitsPath(unit).name(1:end-25),'_F40BG25.bmp']);

end
