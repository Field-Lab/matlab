
cd('/mnt/muench_data/user/alexandra/scripts')


chs=1:60;
chs(15)=[];
keyWord='quick';
date='20121023';
path2data=['/mnt/muench_data/user/alexandra/MEA_data/',date,'/LFP/'];

lfpData=dir([path2data,'*',keyWord,'*.mat']);

hekapath=['/mnt/muench_data/user/alexandra/MEA_data/',date,'/HEKA/'];
heka=dir([hekapath,'*',keyWord,'*.phys']);
% get protocols
whiteFlash=cell(14,1);
blackFlash=cell(14,1);
for i=1:14
    protocol=read_header_field_heka(hekapath, heka(i).name, 'Stimulus Protocol');
    protocol(protocol(:,2)==1,2)=50;
    protocol(protocol(:,2)==0,2)=10;
    protocol=round(protocol(2:end,[1,2]));
    whiteEnds=protocol(:,2)==50;
    whiteBegins=[whiteEnds(2:end); 0];
    whiteFlash{i}=[protocol(logical(whiteBegins),1) protocol(logical(whiteEnds),1)];
    blackEnds=protocol(:,2)==10;
    blackBegins=[blackEnds(2:end); 0];
    blackFlash{i}=[protocol(logical(blackBegins),1) protocol(logical(blackEnds),1)];    
end

load([path2data,lfpData(44).name])
figure
nds='87654321234567'
xl=[-4500 750];
for i=1:12
    subplot(3,4,i)
    whiteLFP=0;
    for j=1:5
        whiteLFP=whiteLFP+lfp(whiteFlash{i}(j,1)-200:whiteFlash{i}(j,1)+3500,i);
    end
    whiteLFP=whiteLFP/5;
    plot(-200:1:3500,whiteLFP,'LineWidth',2)
    axis([-200 3500 xl])
    line([0,0],xl,'color','k')
    line([2000,2000],xl,'color','k')
    title(['ND',nds(i)])
end


% FFFlicker
hekapath=['/mnt/muench_data/user/alexandra/MEA_data/',date,'/HEKA/'];
heka=dir([hekapath,'*.phys']);

file_list=[];
for i=1:length(heka)
    if ~isempty(regexp(heka(i).name,'HCseq', 'once'))
        file_list=[file_list i];
    end
end


protocol=read_header_field_heka(hekapath, heka(file_list(1)).name, 'Stimulus Protocol');
protocol(protocol(:,3)==0&protocol(:,4)==0,4)=-1;
protocol=round(protocol(2:end,[1,4]));
protocol_accum=zeros(size(protocol,1),2,length(file_list));
for i=1:length(file_list)
        protocol=read_header_field_heka(hekapath, heka(file_list(i)).name, 'Stimulus Protocol');
        protocol(protocol(:,3)==0&protocol(:,4)==0,4)=-1;
        protocol=round(protocol(2:end,[1,4]));
        protocol_accum(:,:,i)=protocol;
end

keyWord='HCseq'
path2data=['/mnt/muench_data/user/alexandra/MEA_data/',date,'/LFP/'];
lfpData=dir([path2data,'*',keyWord,'*.mat']);

figure
all_ch=0;
for ch=1:59
    load([path2data,lfpData(ch).name])
    
    HighFilter=zeros(500,4);
    cnt=1;
    for trial_counter=21:24
        %     trial_counter
        flips=protocol_accum(:,1,trial_counter); % flips in ms, rounded
        pxls=((protocol_accum(:,2,trial_counter)-30)/30); % pxs in ms, normalized -1:1
        a=lfp(:,trial_counter);
        
        tmp=zeros(flips(end)+1,1); % 1 column - color, 2 column - spikes
        tmp([1; flips(1:end-1)])=pxls;
        
        % fill stimulus down to ms with brightness values where spikes
        % happen
        for i=1:50
            subscr=(flips(1:end-1)+i);
            subscr(tmp(flips(1:end-1)+i,1)~=0)=[];
            tmp(subscr,1)=tmp(subscr-1,1);
        end
        b=0;
        for i=1:length(a)-2000
            b=b+a(i+500)*tmp(i+499:-1:i);
        end
        b=b/(length(a)-2000);
        HighFilter(:,cnt)=b;
        cnt=cnt+1;        
    end
    all_ch=all_ch+HighFilter;
    subplot(6,10,ch)
    plot(HighFilter,'LineWidth',2)
end
subplot(6,10,60)
plot(all_ch/59,'LineWidth',2)
    



figure
cnt=1;
for i=1:4:48
    subplot(3,4,cnt)
    plot(HighFilter(:,i:i+3),'LineWidth',2)
    title(['ND',nds(cnt)])
    axis([0 500 -70 50])
    cnt=cnt+1;    
end









chs=1:60;
chs(15)=[];
for ch=1:59
    ch
    load([path2data,lfpData(ch).name])
    
    HighFilter=zeros(500,length(file_list));
    for trial_counter=1:length(file_list)
        %     trial_counter
        flips=protocol_accum(:,1,trial_counter); % flips in ms, rounded
        pxls=((protocol_accum(:,2,trial_counter)-30)/30); % pxs in ms, normalized -1:1
        a=lfp(:,trial_counter);
        
        tmp=zeros(flips(end)+1,1); % 1 column - color, 2 column - spikes
        tmp([1; flips(1:end-1)])=pxls;
        
        % fill stimulus down to ms with brightness values where spikes
        % happen
        for i=1:50
            subscr=(flips(1:end-1)+i);
            subscr(tmp(flips(1:end-1)+i,1)~=0)=[];
            tmp(subscr,1)=tmp(subscr-1,1);
        end
        b=0;
        for i=1:length(a)-2000
            b=b+a(i+500)*tmp(i+499:-1:i);
        end
        b=b/(length(a)-2000);
        HighFilter(:,trial_counter)=b;   
    end
    save(['/mnt/muench_data/user/alexandra/MEA_data/',date,'/LFP_filters/lfp_filters_CH',int2str(chs(ch)),'_HCseq1'],'HighFilter')
end





nds='87654321234567'

all_ch=0;
first_last=zeros(59,16);
for ch=1:59
    load(['/mnt/muench_data/user/alexandra/MEA_data/',date,'/LFP_filters/lfp_filters_CH',int2str(chs(ch)),'_HCseq1.mat'])
    tmp=HighFilter(:,sort([1:12:96 12:12:96]));
    [val,ind]=min(tmp);
    first_last(ch,1:16)=ind;
%     all_ch=all_ch+HighFilter;
end

figure
cnt=1;
for i=1:12:96
    subplot(2,4,cnt)
    toPlot=all_ch(:,[i,i+11]);
    toPlot=-toPlot./repmat(min(toPlot),500,1);
    plot(toPlot,'LineWidth',2)
    title(['ND',nds(cnt)])
    axis([0 500 -1.1 1])
    cnt=cnt+1;    
end

figure
cnt=1;
seq='123123123123';
for i=49:60
    subplot(4,3,cnt)
    plot(all_ch(:,i),'LineWidth',2)
    title(['seq',seq(cnt)])
    axis([0 500 -40 20])
    cnt=cnt+1;    
end
clear k p
cnt=1;
for i=1:2:16
    k(cnt,1)=mean(first_last(:,i+1)-first_last(:,i));
    k(cnt,2)=std(first_last(:,i+1)-first_last(:,i))/sqrt(59);
    [~,p(cnt)]=ttest(first_last(:,i+1),first_last(:,i));
    cnt=cnt+1;
end

figure
set(gcf,'position',[1 31 1680 946])
bar(k(:,1))
hold on
errorbar(k(:,1),k(:,2),'.r')
for i=1:8
    text(i-0.3,k(i,1)+k(i,2)+1,[int2str(k(i,1)),'+-',num2str(k(i,2)),'ms']);
    text(i-0.3,k(i,1)+k(i,2)+2.5,['p=',num2str(p(i))]);
end
set(gca,'XTick',1:8,'xticklabel',{'8','7','6','5','4','3','2','1'})
axis([0 9 -30 55])



all_lat=zeros(59,96);
for ch=1:59
    load(['/mnt/muench_data/user/alexandra/MEA_data/',date,'/LFP_filters/lfp_filters_CH',int2str(chs(ch)),'_HCseq1.mat'])
    [val,ind]=min(HighFilter(:,1:96));
    all_lat(ch,:)=ind;
end

k(1:96,1)=mean(all_lat);
k(1:96,2)=std(all_lat)/sqrt(59);
figure
hold on
for i=13:12:96
    plot(i:i+11,k(i:i+11,1),'color','b','linewidth',2)
    errorbar(i:i+11,k(i:i+11,1),k(i:i+11,2),'.','color','r')
end
set(gca,'XTick',18:12:96,'xticklabel',{'7','6','5','4','3','2','1'})
axis([12 97 0 190])



