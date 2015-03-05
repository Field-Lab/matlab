date='20121023'

unitsPath=dir(['/mnt/muench_data/user/alexandra/MEA_data/',date,'/easy_formatted_units/*FFFlash*'])
trials=size(spike_info.name_info,1);
[rows,cols]=opt_subplots(trials*2);
nds='8765432123456';
for unit=1:length(unitsPath)
    load(['/mnt/muench_data/user/alexandra/MEA_data/',date,'/easy_formatted_units/',unitsPath(unit).name])
    white_flash=zeros(trials,4501);
    black_flash=zeros(trials,4501);
    for trial=1:trials
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
    
    for nd=1:trials
        subplot(ceil(trials*2/6),6,cnt)
        plot(white_flash(nd,:)','LineWidth',3);
        line([500,500],[0,max(white_flash(:))],'Color','k','LineWidth',2);
        line([2500,2500],[0,max(white_flash(:))],'Color','k','LineWidth',2);
        axis tight
        set(gca,'XTick',0)
        title(['ND',nds(nd),' Flash RGB 50'],'FontSize',14,'FontWeight','bold','color',[0.6 0.6 0.6])

        subplot(ceil(trials*2/6),6,cnt+1)
        plot(black_flash(nd,:)','LineWidth',3,'color','r');
        line([500,500],[0,max(black_flash(:))],'Color','k','LineWidth',2);
        line([2500,2500],[0,max(black_flash(:))],'Color','k','LineWidth',2);
        axis tight
        set(gca,'XTick',0)
        title(['ND',nds(nd),' Flash RGB 10'],'FontSize',14,'FontWeight','bold')
        cnt=cnt+2;
    end
    subplot('Position',[0.5 0.95 0.00001 0.00001])
    set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
    title([unitsPath(unit).name(1:end-25),'. Spike Rate in response to FF flash 2s, BG30. Each line - mean of 5 flashes'],'FontSize',15,'FontWeight','bold','Interpreter','None')
    saveas(gcf,['F:\',date,'\FFFlash_plots\',unitsPath(unit).name(1:end-25),'_BG30.bmp']);

end
