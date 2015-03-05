date='20121023'

path2data=['S:\user\alexandra\MEA_data\',date,'\'];
unitsPath=dir([path2data,'easy_formatted_units\*FFFlash*']);
path2save=[path2data,'FFFlash_plots\'];

if ~exist(path2save,'file')
    mkdir(path2save);
end
load([path2data,'easy_formatted_units\',unitsPath(1).name])
trials=size(spike_info.name_info,1);
nds='87654321234567';
for unit=1:length(unitsPath)
    load([path2data,'easy_formatted_units\',unitsPath(unit).name])
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
    set(gcf,'position',[1  31 1680 946],'color',[1 1 1]*0.8);
    cnt=1;
    
    for nd=1:trials
        subplot(4,7,cnt)
        plot(white_flash(nd,:)','LineWidth',2);
        hold on
        area(1:500,white_flash(nd,1:500)','FaceColor',[1 1 1]*0.5);
        area(500:2500,white_flash(nd,500:2500)','FaceColor',[1 1 1]*0.8);
        area(2500:4501,white_flash(nd,2500:4501)','FaceColor',[1 1 1]*0.5);
        line([500,500],[0,max(white_flash(:))],'Color','b','LineWidth',1);
        line([2500,2500],[0,max(white_flash(:))],'Color','b','LineWidth',1);
        axis tight
        set(gca,'XTick',0)
        title(['ND',nds(nd),' On Flash (50)'],'FontSize',12,'Backgroundcolor',[1 1 1]*0.8)
        cnt=cnt+1;
    end
    for nd=1:trials
        subplot(4,7,cnt)
        plot(black_flash(nd,:)','LineWidth',2);
        hold on
        area(1:500,black_flash(nd,1:500)','FaceColor',[1 1 1]*0.5);
        area(500:2500,black_flash(nd,500:2500)','FaceColor',[1 1 1]*0.2);
        area(2500:4501,black_flash(nd,2500:4501)','FaceColor',[1 1 1]*0.5);
        line([500,500],[0,max(black_flash(:))],'Color','b','LineWidth',1);
        line([2500,2500],[0,max(black_flash(:))],'Color','b','LineWidth',1);
        axis tight
        set(gca,'XTick',0)
        title(['ND',nds(nd),' Off Flash (10)'],'FontSize',12,'Backgroundcolor',[1 1 1]*0.5)
        cnt=cnt+1;
    end
    subplot('Position',[0.5 0.95 0.00001 0.00001])
    set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
    title([unitsPath(unit).name(1:end-25),'. Spike Rate, Full Field Flash 2s, BG 30, end of ND. Mean of 5 flashes'],'FontSize',15,'FontWeight','bold','Interpreter','None')
    saveas(gcf,[path2save,unitsPath(unit).name(1:end-25),'_BG30.emf']);

end
