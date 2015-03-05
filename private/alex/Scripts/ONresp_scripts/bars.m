clear
date='20130703a'
% date='20130703b'

codeWord='bars';

mainpath=['/Users/alexth/Desktop/old_stuff/',date,'/'];
hekapath=['/Users/alexth/Desktop/old_stuff/',date,'/HEKA/'];
heka=dir([hekapath,'*.phys']);

path2save=['/Users/alexth/Desktop/old_stuff/ONresp/',date,'/'];
if ~exist(path2save,'dir')
    mkdir(path2save);
end

% select heka files
clear ndFullList
file_list=[];
cnt=1;
for i=1:length(heka)
    if ~isempty(regexp(heka(i).name,codeWord, 'once'))
        file_list=[file_list i];
        
        a=regexp(heka(i).name,'ND');
        a=heka(i).name(a+2);
        if length(a)>1
            ndFullList(cnt)=str2num(a(1))+str2num(a(2));
        else
            ndFullList(cnt)=str2num(a(1));
        end
        cnt=cnt+1;
    end
end


% get protocol times
protocol=read_header_field_heka(hekapath, heka(file_list(1)).name, 'Stimulus Protocol');
protocol=round(protocol(2:end-1,[1,2,3]));

units=dir([mainpath,'units/*.mat']);
bars=zeros(4500,length(file_list)*8,length(units));

names=cell(length(units),1);
barcode=zeros(length(file_list)*8,2);

for cnt=1:length(units)
    
    load([mainpath,'units/',units(cnt).name]);
    trialNumber=1;
    
    for i=1:length(file_list)
        protocol=read_header_field_heka(hekapath, heka(file_list(i)).name, 'Stimulus Protocol');
        protocol=round(protocol(2:end-1,[1,2,3]));
        spikes=round(cell2mat(unit{1,2}(file_list(i),2))*1000/cell2mat(unit{1,1}(5,2))); %now in ms
        conv=convolved(spikes,40,102000);
        conv=conv(121:end-120);
        for j=1:2:size(protocol,1)-1
            bars(:,trialNumber,cnt)=conv(protocol(j,1)-499:protocol(j,1)+4000);
            barcode(trialNumber,1:2)=protocol(j+1,2:3);
            trialNumber=trialNumber+1;
        end
        
    end  
    
    if isempty(regexp(units(cnt).name,date, 'once'))
        names{cnt}=[units(cnt).name(1),date,units(cnt).name(10:end-4)];
    else
        names{cnt}=units(cnt).name(1:end-4);
    end
end
barcode(barcode==2)=-1;
barcode(barcode==3)=1;
save([path2save,'bars'],'ndFullList','bars','names','barcode')



%% Plot

clear

date='20130703a'
path2load=['/Users/alexth/Desktop/old_stuff/ONresp/',date,'/'];
load([path2load,'bars'])
load([path2load,'quick'])


coord_pos=zeros(8,10);
coord_neg=zeros(8,10);
coord=-1000:200:400;
for i=1:8
    a=barcode(:,2)==coord(i);
    coord_pos(i,1:10)=find(a&barcode(:,1)==1);
    coord_neg(i,1:10)=find(a&barcode(:,1)==-1);
end

path2save=['/Users/alexth/Desktop/old_stuff/ONresp/',date,'/bars_plot/'];
if ~exist(path2save,'dir')
    mkdir(path2save);
end

wid=0.10;

figure
set(gcf,'position',[2         399        1283         538])
for k=1:size(bars,3)

    for i=1:8
        subplot('position',[0.015+(i-1)*wid*1.1 0.52 wid 0.4])
        a=mean(bars(1:4000,coord_neg(i,1:5),k),2);
        b=mean(bars(1:4000,coord_pos(i,1:5),k),2);
        hold off
        plot(a)
        hold on
        plot(b,'r')
        title(['ND6, ',int2str(coord(i))])
        
        line([500 500],[0 80],'color','k')
        line([2500 2500],[0 80],'color','k')
        axis([1 4000 0 80])
        set(gca,'xtick',0,'xticklabel','','yticklabel','')
        
        subplot('position',[0.015+(i-1)*wid*1.1 0.02 wid 0.4])
        a=mean(bars(1:4500,coord_neg(i,6:10),k),2);
        b=mean(bars(1:4500,coord_pos(i,6:10),k),2);
        hold off
        plot(a,'g')
        hold on
        plot(b,'m')
        title(['ND4, ',int2str(coord(i))])
        line([500 500],[0 80],'color','k')
        line([2500 2500],[0 80],'color','k')
        axis([1 4000 0 80])
        set(gca,'xtick',0,'xticklabel','','yticklabel','')
    end
    
    i=9;
    subplot('position',[0.015+(i-1)*wid*1.1 0.52 wid 0.4])
    hold off
    plot(mean(black_flash(:,17:26,k),2),'b')
    hold on
    plot(mean(white_flash(:,17:26,k),2),'r')
    title('ND6 full field')
    line([500 500],[0 80],'color','k')
    line([2500 2500],[0 80],'color','k')
    axis([1 4000 0 80])
    set(gca,'xtick',0,'xticklabel','','yticklabel','')
            
    subplot('position',[0.015+(i-1)*wid*1.1 0.02 wid 0.4])
    hold off
    plot(mean(black_flash(:,35:44,k),2),'g')
    hold on
    plot(mean(white_flash(:,35:44,k),2),'m')
    title('ND4 full field')
    line([500 500],[0 80],'color','k')
    line([2500 2500],[0 80],'color','k')
    axis([1 4000 0 80])
    set(gca,'xtick',0,'xticklabel','','yticklabel','')
        
    
    subplot('Position',[0.5 0.95 0.00001 0.00001])
    set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
    title([names{k},'  BG 10, BG50'],'FontSize',12,'FontWeight','bold','Interpreter','None')
    saveas(gcf,[path2save,names{k},'.png'])
end



%%%%% ND6 %%%%%


path2save=['/Users/alexth/Desktop/old_stuff/ONresp/',date,'/bars_plot_for_Thomas/'];
if ~exist(path2save,'dir')
    mkdir(path2save);
end


figure
set(gcf,'position',[40         139        1202         788])
for k=1:size(bars,3)

    for i=1:8
        subplot(4,2,i)
        a=mean(bars(:,coord_pos(i,1:5),k),2);
        b=mean(bars(:,coord_neg(i,1:5),k),2);
        hold off
        plot(a,'m','linewidth',2)
        hold on
        plot(4701:9200, b,'c','linewidth',2)
        plot(mean(white_flash(:,17:26,k),2),'r')
        plot(4701:9200, mean(black_flash(:,17:26,k),2),'b')
        legend('white bar','black bar','white FF','black FF','location','bestoutside')
         title(['ND6, position: ',int2str(coord(i))])
        
        axis tight
        line([500 500],[0 80],'color','k')
        line([2500 2500],[0 80],'color','k')  
        line([500 500]+4700,[0 80],'color','k')
        line([2500 2500]+4700,[0 80],'color','k')  
    end

    subplot('Position',[0.5 0.95 0.00001 0.00001])
    set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
    title([names{k},'  ND6, BG 50, BG 10'],'FontSize',12,'FontWeight','bold','Interpreter','None')
    saveas(gcf,[path2save,names{k},'_ND6.png'])
end