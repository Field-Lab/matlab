clear
date='20130226'

        
path2data=['S:\data\alexandra\MEA_data\',date,'\easy_formatted_units\'];
flashes=dir([path2data,'*FFFlash*']);
moments=dir([path2data,'*Moment*']);
nds=dir([path2data,'*nd*']);


path2save=['S:\data\alexandra\MEA_data\',date,'\full_info\'];

if ~exist(path2save,'file')
    mkdir(path2save);
end
load([path2data,flashes(1).name])
trials=size(spike_info.name_info,1);


filterspath=['S:\data\alexandra\MEA_data\',date,'\FFFlicker_LF\'];
filters=dir([filterspath,'*FFFlicker_linear_filter.mat']);


for unit=1:length(flashes)
    
    close(figure(1))
    figure(1)
    set(gcf,'position',[1  31 1680 946],'color',[1 1 1]*0.8);  
    
    %nd
    load([path2data,nds(unit).name])
    
    kk=1;
    for i=1:9
        subplot(9,4,kk)
        spikes=cell2mat(spike_info.spike_times(10-i));
        flips=cell2mat(spike_info.flip_times(10-i));
        conv=convolved(spikes,40,flips(end,1));
        conv=conv(121:end-120);        
        plot(conv,'Linewidth',2)
        line([flips(1,1)-1000,flips(1,1)-1000],[0,max(conv)],'color','k')
        line([flips(6,1),flips(6,1)],[0,max(conv)],'color','k')
        title(['ND',int2str(i)])        
        axis tight
%         set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
        kk=kk+4;
    end
    
    
    % FULL FIELD FLASH
    load([path2data,flashes(unit).name])
    load([filterspath,filters(unit).name])
    LF=HighFilters';
    
    clear white_flash black_flash
    white_flash_tmp=zeros(trials,4501);
    black_flash_tmp=zeros(trials,4501);
    for trial=1:trials
        spikes=cell2mat(spike_info.spike_times(trial));
        flips=cell2mat(spike_info.flip_times(trial));
        flips=flips(:,1);
        conv=convolved(spikes,40,flips(end));
        conv=conv(121:end-120);
        for st=1:4:20
            white_flash_tmp(trial,1:4501)=white_flash_tmp(trial,1:4501)+conv(flips(st)-500:flips(st)+4000);
            black_flash_tmp(trial,1:4501)=black_flash_tmp(trial,:)+conv(flips(st+2)-500:flips(st+2)+4000);
        end
    end
    
   
    cnt=1;
    for i=1:9:81
        white_flash(cnt,1:4501)=mean(white_flash_tmp(i:i+8,:));
        black_flash(cnt,1:4501)=mean(black_flash_tmp(i:i+8,:));
        cnt=cnt+1;
    end
    white_flash=white_flash/5;
    black_flash=black_flash/5;
    
    maxx=max(max([white_flash(:,1:4500) black_flash(:,1:4500)]));
    kk=2;
    for i=1:9
        subplot(9,4,kk)        
        hold on       
        area(1:500,white_flash(10-i,1:500)','FaceColor',[1 1 1]*0.5);
        area(500:2500,white_flash(10-i,500:2500)','FaceColor',[1 1 1]*0.8);
        area(2500:4500,white_flash(10-i,2500:4500)','FaceColor',[1 1 1]*0.5);
       
        area(4501:5000,black_flash(10-i,1:500)','FaceColor',[1 1 1]*0.5);
        area(5001:7000,black_flash(10-i,501:2500)','FaceColor',[1 1 1]*0.1);
        area(7001:9000,black_flash(10-i,2501:4500)','FaceColor',[1 1 1]*0.5);
        plot([white_flash(10-i,1:4500) black_flash(10-i,1:4500)],'r','LineWidth',2);
        
        line([500,500],[0,maxx],'Color','k','LineWidth',2);
        line([2500,2500],[0,maxx],'Color','k','LineWidth',2);
        line([4500,4500],[0,maxx],'color','k','LineWidth',3);
        line([5001,5001],[0,maxx],'Color','k','LineWidth',2);
        line([7001,7001],[0,maxx],'Color','k','LineWidth',2);
        axis tight
        set(gca,'XTick',0,'XTickLabel','')
        title(['ND',int2str(i)])
        kk=kk+4;
    end
    
    
    % Moment    
    load([path2data,moments(unit).name])
    clear white_flash black_flash
    white_flash_tmp=zeros(trials,2501);
    black_flash_tmp=zeros(trials,2501);
    for trial=1:trials
        spikes=cell2mat(spike_info.spike_times(trial));
        flips=cell2mat(spike_info.flip_times(trial));
        flips=flips(:,1);
        conv=convolved(spikes,40,flips(end));
        conv=conv(121:end-120);
        for st=1:4:20
            white_flash_tmp(trial,1:2501)=white_flash_tmp(trial,:)+conv(flips(st)-500:flips(st)+2000);
            black_flash_tmp(trial,1:2501)=black_flash_tmp(trial,:)+conv(flips(st+2)-500:flips(st+2)+2000);
        end
    end
    
   
    cnt=1;
    for i=1:9:81
        white_flash(cnt,1:2501)=mean(white_flash_tmp(i:i+8,:));
        black_flash(cnt,1:2501)=mean(black_flash_tmp(i:i+8,:));
        cnt=cnt+1;
    end
    white_flash=white_flash/5;
    black_flash=black_flash/5;
   
    maxx=max(max([white_flash(:,1:2500) black_flash(:,1:2500)]));
    kk=3;
    for i=1:9
        subplot(9,4,kk)        
        hold on       
        plot([white_flash(10-i,1:2500) black_flash(10-i,1:2500)],'r','LineWidth',2);

        axis tight
        set(gca,'XTick',0,'XTickLabel','')
        title(['ND',int2str(i)])
        kk=kk+4;
    end

    % LF
    kk=4;
    for i=1:9
        subplot(9,4,kk)
        plot(LF(9*(10-i-1)+1:9*(10-i),:)')
        line([0 500],[0 0],'color',[1 1 1]*0.5)
        title(['ND',int2str(i)])        
        axis tight
        set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
        kk=kk+4;
    end
    
    
    
    subplot('Position',[0.5 0.97 0.00001 0.00001])
    set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
    title(flashes(unit).name(1:end-4),'FontSize',12,'FontWeight','bold','Interpreter','None')
    saveas(gcf,[path2save,flashes(unit).name(1:end-25),'.emf']);

end
        
 

