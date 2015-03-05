clear
date='20121004'

unitsPath=dir(['S:\data\alexandra\MEA_data\',date,'\easy_formatted_units\*chirp*'])
load(['S:\data\alexandra\MEA_data\',date,'\easy_formatted_units\',unitsPath(1).name])
trials=size(spike_info.name_info,1);

tmp3=spike_info.flip_times{1}(2:end-1,2);
[a,b]=findpeaks(tmp3);
c=spike_info.flip_times{1}(b+1,1);
c=c/10000;

nds='7654432123456';
close(figure(1))
figure(1)
set(gcf,'position',[1  31 1680 946]);
for unit=1:length(unitsPath)
    load(['S:\data\alexandra\MEA_data\',date,'\easy_formatted_units\',unitsPath(unit).name])
    chirp=zeros(trials,25000);
    for trial=1:trials
        spikes=cell2mat(spike_info.spike_times(trial));
        flips=cell2mat(spike_info.flip_times(trial));
        flips=flips(:,1);
        conv=convolved(spikes,40,flips(end,1));
        conv=conv(121:end-120);
        chirp(trial,1:length(conv))=conv;
    end

    cnt=1;
    
    for nd=11:5:45
        subplot(7,1,cnt)
        hold off
        tmp=mean(chirp(nd:nd+4,1:20000));
        plot((1:length(tmp))/10000,tmp,'LineWidth',1);
        k=mean(tmp(500:1500));
        line([0,length(tmp)/10000],[k,k],'color','m')
        title(['ND',nds(cnt)],'FontSize',12,'FontWeight','bold')
        for i=1:length(c)
            line([c(i),c(i)],[0,max(tmp)],'color','k')
        end        
        hold on
        plot((1:length(tmp))/10000,tmp,'LineWidth',2);
        axis tight
        set(gca,'xtick',0)
        cnt=cnt+1;
    end
    subplot('Position',[0.5 0.95 0.00001 0.00001])
    set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
    title([unitsPath(unit).name(1:end-23),'. Chirp spike rate, average of 5 trials, end of ND'],'FontSize',16,'FontWeight','bold','Interpreter','None')
    saveas(gcf,['S:\data\alexandra\MEA_data\',date,'\chirp_plots\',unitsPath(unit).name(1:end-23),'_chirp.emf']);

end
