clear
date='20120928'

path2fit=['S:\user\alexandra\MEA_data\',date,'\fit_res_HC\'];
fits=dir([path2fit,'\*.mat']);

        
path2data=['S:\user\alexandra\MEA_data\',date,'\'];
unitsPath=dir([path2data,'easy_formatted_units\*FFFlash*']);
path2save=[path2data,'Full_info\'];

if ~exist(path2save,'file')
    mkdir(path2save);
end
load([path2data,'easy_formatted_units\',unitsPath(1).name])
% trials=size(spike_info.name_info,1);
trials=7;
% nds='876543212345';
nds='7654321';


chirpPath=dir(['S:\user\alexandra\MEA_data\',date,'\easy_formatted_units\*chirp*']);
load(['S:\user\alexandra\MEA_data\',date,'\easy_formatted_units\',chirpPath(1).name])
trialsChirp=size(spike_info.name_info,1);
flips=spike_info.flip_times{1}(:,1); % flips in ms, rounded
pxls=spike_info.flip_times{1}(:,2); % pxs in ms, normalized -1:1
pxls(pxls==-1)=30;
tmp1=zeros(flips(end)+1,1); % 1 column - color, 2 column - spikes
tmp1([1; flips(1:end-1)],1)=pxls;
    
for i=1:50
    subscr=(flips(1:end-1)+i);
    subscr(tmp1(flips(1:end-1)+i,1)~=0)=[];
    tmp1(subscr,1)=tmp1(subscr-1,1);
end
tmp1(tmp1==0)=30;
cdata = tmp1/60*0.8+0.1;
cdata=repmat(cdata',400,1);

filterspath=['F:\',date,'\FFFlicker_LF\'];
filters=dir([filterspath,'*FFFlicker_linear_filter.mat']);

for unit=1:length(fits)
    load([path2fit,fits(unit).name]);
    name=str2num(fits(unit).name(end-18:end-15));
    if name==22
        break
    end
end

for unit=1:length(fits)
    load([path2fit,fits(unit).name]);
    name=fits(unit).name(1:end-15);
    cellType=sum(sign(common_res_fit(1,:)));
    
    % FULL FIELD FLASH
    load([path2data,'easy_formatted_units\',name,'___spike_info_FFFlash'])
    
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
    white_flash=white_flash(1:7,:)/5;
    black_flash=black_flash(1:7,:)/5;
    close(figure(1))
    figure(1)
    set(gcf,'position',[1  31 1680 946],'color',[1 1 1]*0.8);
    
     maxx=max(max([white_flash(:,1:4500) black_flash(:,1:4500)]));
    for i=1:7
        subplot('Position',[0.7 0.01+(0.14*(i-1)) 0.15 0.12])        
        hold on
       
        area(1:500,white_flash(i,1:500)','FaceColor',[1 1 1]*0.5);
        area(500:2500,white_flash(i,500:2500)','FaceColor',[1 1 1]*0.8);
        area(2500:4500,white_flash(i,2500:4500)','FaceColor',[1 1 1]*0.5);
       
        area(4501:5000,black_flash(i,1:500)','FaceColor',[1 1 1]*0.5);
        area(5001:7000,black_flash(i,501:2500)','FaceColor',[1 1 1]*0.1);
        area(7001:9000,black_flash(i,2501:4500)','FaceColor',[1 1 1]*0.5);
        plot([white_flash(i,1:4500) black_flash(i,1:4500)],'r','LineWidth',2);
        
        line([500,500],[0,maxx],'Color','k','LineWidth',2);
        line([2500,2500],[0,maxx],'Color','k','LineWidth',2);
        line([4500,4500],[0,maxx],'color','k','LineWidth',3);
        line([5001,5001],[0,maxx],'Color','k','LineWidth',2);
        line([7001,7001],[0,maxx],'Color','k','LineWidth',2);
        axis tight
    end

    % CHIRP
    
    load(['S:\user\alexandra\MEA_data\',date,'\easy_formatted_units\',name,'___spike_info_chirp'])
    chirp=zeros(trialsChirp,25000);
    for trial=1:trialsChirp
        spikes=cell2mat(spike_info.spike_times(trial));
        flips=cell2mat(spike_info.flip_times(trial));
        flips=flips(:,1);
        conv=convolved(spikes,40,flips(end,1));
        conv=conv(121:end-120);
        chirp(trial,1:length(conv))=conv;        
    end
    
    cnt=1;
    for nd=1:5:35
        tmp(cnt,1:20000)=mean(chirp(nd:nd+4,1:20000));
        cnt=cnt+1;
    end
    
    for i=1:7
        subplot('Position',[0.03 0.01+(0.14*(i-1)) 0.65 0.12])
%         imagesc(cdata);
        set(gca,'YDir','normal')
%         colormap(gray);
        hold on
        for nd=(i-1)*5+1:(i-1)*5+5
            plot(1:20000,chirp(nd,1:20000),'b','LineWidth',2);
        end
        plot(1:20000,tmp(i,:),'r','LineWidth',2);
        line([1,20000],[mean(tmp(i,500:1500)) mean(tmp(i,500:1500))],'color','k')
        axis([2000 20000 0 max(tmp(:))]);
        ylabel(['ND',nds(i)])
    end
    
    
    % FLICKER
    
    load([filterspath,name,'_FFFlicker_linear_filter'])
%     HighFilters=HighFilters(:,13:end);
%     HighSpikeCount=HighSpikeCount(13:end);
    minn=min(min(HighFilters(:,1:2:end)));
    maxx=max(max(HighFilters(:,1:2:end)));
    i=1;
    for cnt=1:12:12*7
        subplot('Position',[0.87 0.01+(0.14*(i-1)) 0.12 0.12])
        goto=min(cnt+11,size(HighFilters,2));
        plot(HighFilters(:,cnt:2:goto))
        line([0 500],[0 0],'color',[1 1 1]*0.5)
        hleg=legend(int2str(HighSpikeCount(cnt:2:goto)'));
        legend('boxoff')
        set(hleg,'Fontsize',6)
        axis([0 500 minn maxx])
        i=i+1;
    end
    
    
    
    % COMMON TITLE
    if cellType>0 %ON cell
        lab='ON';
    else
        lab='OFF';
    end
    subplot('Position',[0.5 0.97 0.00001 0.00001])
    set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
    title([name, '     ', lab, ' cell'],'FontSize',12,'FontWeight','bold','Interpreter','None')
    saveas(gcf,[path2save,name,'.emf']);


end
    