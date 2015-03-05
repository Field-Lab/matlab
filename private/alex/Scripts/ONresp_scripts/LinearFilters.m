cd('/Users/alexth/Desktop/scripts/ONresp_scripts')

%% Prepare matrix for Full Field Flicker
clear
date='20131025a'
codeWord='HC'

mainpath=['/Users/alexth/Desktop/old_stuff/',date,'/'];
hekapath=['/Users/alexth/Desktop/old_stuff/',date,'/HEKA/'];
heka=dir([hekapath,'*.phys']);

path2save=['/Users/alexth/Desktop/old_stuff/ONresp/',date,'/'];
if ~exist(path2save,'dir')
    mkdir(path2save);
end

% select heka files
file_list=[];
for i=1:length(heka)
    if ~isempty(regexp(heka(i).name,codeWord, 'once'))
        file_list=[file_list i];
    end
end

% get protocol

protocol=read_header_field_heka(hekapath, heka(file_list(2)).name, 'Stimulus Protocol');
protocol(protocol(:,3)==0&protocol(:,4)==0,4)=-1;
protocol=round(protocol(2:end-2,[1,4]));
correctedProtocols=zeros(size(protocol,1),2,length(file_list));
tmp=(0:size(correctedProtocols,1)-1)'*(7.85/3600);
ndFullList=zeros(1,length(file_list));

for i=1:length(file_list)
    hekaFile=heka(file_list(i)).name;
    protocol=read_header_field_heka(hekapath, hekaFile, 'Stimulus Protocol');
    protocol(protocol(:,3)==0&protocol(:,4)==0,4)=-1;
    protocol=round(protocol(2:end-2,[1,4]));
    protocol(:,1)=protocol(:,1)+tmp;
    correctedProtocols(:,:,i)=round(protocol);
    
    a=regexp(hekaFile,'ND');
    a=hekaFile(a+2);
    a=a(isstrprop(a, 'digit'));
    for t=1:length(a)
        ndFullList(i)=ndFullList(i)+str2num(a(t));
    end
end

clear protocol

startTimes=reshape(correctedProtocols(1,1,:),1,length(file_list));
endTimes=reshape(correctedProtocols(end,1,:),1,length(file_list));
maxLength=max(endTimes-startTimes);
flicker=zeros(maxLength,length(file_list));

for i=1:length(file_list)
%     if i==146
%         correctedProtocols(:,1,i)=correctedProtocols(:,1,i-1);
%     end
    flips=correctedProtocols(:,1,i);
    pxls=(correctedProtocols(2:end,2,i)-30)/30+256;
    
    flips=flips-flips(1)+1;
    
    tmp=zeros(flips(end),1);
    tmp(flips(2:end),1)=pxls;
    
    % fill stimulus down to ms with brightness values where spikes
    % happen
    droppedMax=max(diff(flips));
    for j=1:droppedMax+1
        subFrame=flips(2:end)-j;
        subFrame(subFrame<1)=[];
        subFrame(find(tmp(subFrame,1)))=[];
        tmp(subFrame,1)=tmp(subFrame+1,1);
    end
    tmp=tmp-256;
    flicker(1:length(tmp),i)=tmp;
end

save([path2save,'protocols_',codeWord],'ndFullList','mainpath','date','flicker','startTimes','maxLength','correctedProtocols','codeWord','file_list')



%% Get filter

clear
date='20131025a'
codeWord='HC'
path2save=['/Users/alexth/Desktop/old_stuff/ONresp/',date,'/'];
load([path2save,'protocols_',codeWord])

per=1800; % frames to take for calculations
filter_length=500;

units=dir([mainpath, 'units/*.mat']);

names=cell(1,length(units));
LinearFilter=zeros(filter_length,length(file_list),length(units));
SpikeCount=zeros(length(file_list),length(units));
if length(unique(startTimes))>1
    for i=1:size(correctedProtocols,3)
        correctedProtocols(:,1,i)=correctedProtocols(:,1,i)-startTimes(i);
    end
else
    correctedProtocols(:,1,:)=correctedProtocols(:,1,:)-startTimes(1);
end


for cnt=1:length(units)
    
    load([mainpath,'units/',units(cnt).name]);
    for i=1:length(file_list)
        spikes=round(cell2mat(unit{1,2}(file_list(i),2))*1000/cell2mat(unit{1,1}(5,2))); %now in ms
        if ~isempty(spikes)
            spikes=spikes-startTimes(i)+1;
            startPoint=correctedProtocols(1,1,i)+1;
            endPoint=correctedProtocols(per,1,i);
            spikesTMP=spikes(spikes>startPoint+500&spikes<endPoint);
            
            n=zeros(length(spikesTMP),filter_length);
            for k=1:filter_length
                n(:,k)=flicker(spikesTMP-k+1,i);
            end
            LinearFilter(1:filter_length,i,cnt)=sum(n);
            SpikeCount(i,cnt)=length(spikesTMP);
        end
        
        if isempty(regexp(units(cnt).name,date, 'once'))
            names{cnt}=[units(cnt).name(1),date,units(cnt).name(10:end-4)];
        else
            names{cnt}=units(cnt).name(1:end-4);
        end
    end
end

save([path2save,'LF_',codeWord],'ndFullList','LinearFilter','SpikeCount','names','units','per')


%% Plot Linear Filters RInger Switch
clear
date='20130906'
path2load=['/Users/alexth/Desktop/old_stuff/ONresp/',date,'/'];
load([path2load,'LF_HC'])

path2save=[path2load,'/LF_plot_HC/'];
if ~exist(path2save,'dir')
    mkdir(path2save);
end


figure(1)
set(gcf,'Position', [1 352 1011 625])
nds='87654';
for un=1:size(LinearFilter,3)
    for cnt=1:5
        subplot(2,3,cnt)
        plot(LinearFilter(:,ndFullList==str2num(nds(cnt)),un))
        title(['ND',nds(cnt)])
    end
    for i=1:5
        subplot(2,3,i)
        line([0 500],[0 0],'color','k')
        k=LinearFilter(:,:,un);
        axis([0 500 min(k(:)) max(k(:))])
    end
    subplot('Position',[0.5 0.95 0.00001 0.00001])
    set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
    title([names{un},'  HIGH contrast'],'FontSize',12,'FontWeight','bold','Interpreter','None')
    saveas(gcf,[path2save,names{un},'.png'])
    
end


% alternative plot
figure(1)
set(gcf,'Position', [1 352 1011 625])
nds='87654';
for un=1:size(LinearFilter,3)

    subplot(2,3,1)
    plot(LinearFilter(:,1:8,un),'b')
    title('ND8')
    
    subplot(2,3,2)
    hold off
    plot(LinearFilter(:,9:12,un),'b')
    hold on
    plot(LinearFilter(:,13:16,un),'r')
    title('ND7')
    
    subplot(2,3,3)
    hold off
    plot(LinearFilter(:,17:22,un),'b')
    hold on
    plot(LinearFilter(:,23:26,un),'r')
    title('ND6')
    
    subplot(2,3,4)
    hold off
    plot(LinearFilter(:,27:30,un),'b')
    hold on
    plot(LinearFilter(:,31:34,un),'r')
    title('ND5')
    
    subplot(2,3,5)
    hold off
    plot(LinearFilter(:,35:40,un),'b')
    hold on
    plot(LinearFilter(:,41:44,un),'r')
    title('ND4')
    
    for i=1:5
        subplot(2,3,i)
        line([0 500],[0 0],'color','k')
        k=LinearFilter(:,:,un);
        axis([0 500 min(k(:)) max(k(:))])
    end
    subplot('Position',[0.5 0.95 0.00001 0.00001])
    set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
    title([names{un},'  HIGH contrast'],'FontSize',12,'FontWeight','bold','Interpreter','None')
    saveas(gcf,[path2save,names{un},'.png'])
    
end


%% Firing rate

clear

date='20130703b'
codeWord='HC'
path2save=['/Users/alexth/Desktop/old_stuff/ONresp/',date,'/'];
load([path2save,'protocols_',codeWord])
load([path2save,'LF_',codeWord])
FiringRate=zeros(32000,size(SpikeCount,1),size(SpikeCount,2));
for cnt=1:length(units)
    cnt
    load([mainpath,'units/',units(cnt).name]);
    for i=1:length(file_list)
        spikes=round(cell2mat(unit{1,2}(file_list(i),2))*1000/cell2mat(unit{1,1}(5,2))); %now in ms
        if ~isempty(spikes)
            tmp=convolved(spikes,40,32000);
            FiringRate(:,i,cnt)=tmp(121:end-120);
        end
    end
end

nds='87654';
for cnt=1:size(FiringRate,3)
    for i=1:5
        kk=ndFullList==str2num(nds(i));
        tmp=FiringRate(:,kk,cnt);
        tmp=reshape(tmp(2501:31000,:),28500*sum(kk),1);
        
        frSTD_HC(i,cnt)=std(tmp);
        frMean_HC(i,cnt)=mean(tmp);
        
        a=FiringRate(200:1950,kk,cnt);
        frspont(i,cnt)=mean(a(:));
        frSTDspont(i,cnt)=std(a(:));
    end
end


% 20130703a
onOff=[0 1 -1 -1 1 -1 1 1 1 0 1 1 1 1 0 1 0 1 -1 1 -1 1 1];

% 20130703b
onOff=[1 0 1 1 1 1 0 0 1 0 1 -1 0 0 -1 -1 0 1 1 -1 0 -1 1 1 0 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0];


figure
nds='87654';
for i=1:5
    subplot(2,3,i)
    plot(frspont(i,onOff>0),frMean_HC(i,onOff>0)-frspont(i,onOff>0),'*r')
    hold on
    plot(frspont(i,onOff<0),frMean_HC(i,onOff<0)-frspont(i,onOff<0),'*b')
    title(['ND',nds(i)])
    line([0 45],[0 0],'color','k')
end









%% Plot Linear Filters APB
clear
date='2013-09-03'
path2load=['/Users/alexth/Desktop/old_stuff/ONresp/',date,'/'];
load([path2load,'LF_HC'])

path2save=[path2load,'/LF_plot_HC/'];
if ~exist(path2save,'dir')
    mkdir(path2save);
end


figure(1)
set(gcf,'Position', [1 352 1011 625])
nds='7654';
for un=1:size(LinearFilter,3)

    subplot(2,2,1)
    plot(mean(LinearFilter(:,3:6,un),2),'b','linewidth',2)
    title('ND7')
    
    subplot(2,2,2)
    hold off
    plot(mean(LinearFilter(:,7:10,un),2),'b','linewidth',2)
    hold on
    plot(mean(LinearFilter(:,15:18,un),2),'r','linewidth',2)
    plot(mean(LinearFilter(:,27:30,un),2),'g','linewidth',2)
    title('ND6')
    
    subplot(2,2,3)
    hold off    
    plot(mean(LinearFilter(:,31:34,un),2),'b','linewidth',2)
    hold on
    plot(mean(LinearFilter(:,39:42,un),2),'r','linewidth',2)
    plot(mean(LinearFilter(:,51:54,un),2),'g','linewidth',2)
    title('ND5')
    
    subplot(2,2,4)
    hold off
    plot(mean(LinearFilter(:,55:58,un),2),'b','linewidth',2)
    hold on
    plot(mean(LinearFilter(:,63:66,un),2),'r','linewidth',2)
    plot(mean(LinearFilter(:,75:78,un),2),'g','linewidth',2)
    title('ND4')
    
    
%     for i=1:4
%         subplot(2,2,i)
%         line([0 500],[0 0],'color','k')
%         k=LinearFilter(:,:,un);
%         axis([0 500 min(k(:)) max(k(:))])
%     end
    
    subplot('Position',[0.5 0.95 0.00001 0.00001])
    set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
    title([names{un},'  HIGH contrast'],'FontSize',12,'FontWeight','bold','Interpreter','None')
    saveas(gcf,[path2save,names{un},'.png'])
    
end



% special for 20130723a

figure(1)
set(gcf,'Position', [1 352 1011 625])
nds='7654';
for un=1:size(LinearFilter,3)

    subplot(2,2,1)
    plot(mean(LinearFilter(:,3:6,un),2),'b','linewidth',2)
    title('ND7')
    
    subplot(2,2,2)
    hold off
    plot(mean(LinearFilter(:,7:10,un),2),'b','linewidth',2)
    hold on
    plot(mean(LinearFilter(:,15:18,un),2),'r','linewidth',2)
    plot(mean(LinearFilter(:,24:27,un),2),'g','linewidth',2)
    title('ND6')
    
    subplot(2,2,3)
    hold off    
    plot(mean(LinearFilter(:,28:30,un),2),'b','linewidth',2)
    hold on
    plot(mean(LinearFilter(:,35:38,un),2),'r','linewidth',2)
    plot(mean(LinearFilter(:,47:50,un),2),'g','linewidth',2)
    title('ND5')
    
    subplot(2,2,4)
    hold off
    plot(mean(LinearFilter(:,51:54,un),2),'b','linewidth',2)
    hold on
    plot(mean(LinearFilter(:,59:62,un),2),'r','linewidth',2)
    plot(mean(LinearFilter(:,71:74,un),2),'g','linewidth',2)
    title('ND4')
    
    
%     for i=1:4
%         subplot(2,2,i)
%         line([0 500],[0 0],'color','k')
%         k=LinearFilter(:,:,un);
%         axis([0 500 min(k(:)) max(k(:))])
%     end
    
    subplot('Position',[0.5 0.95 0.00001 0.00001])
    set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
    title([names{un},'  HIGH contrast'],'FontSize',12,'FontWeight','bold','Interpreter','None')
    saveas(gcf,[path2save,names{un},'.png'])
    
end








figure(3)
set(gcf,'Position', [1 352 1011 625])
for un=1:33
    
    
    tmp(1:500,1)=mean(LinearFilter(:,(7:10),un),2);
    
    tmp(1:500,2)=mean(LinearFilter(:,(24:27),un),2);
    
    subplot(6,6,un)
    plot(tmp)
    axis([0 500 -100 100])
    
end




figure(4)
set(gcf,'Position', [1 352 1011 625])
for un=1:33
    
    
    tmp(1:500,1)=mean(LinearFilter(:,(1:6),un),2);
    
    
    subplot(6,6,un)
    plot(tmp)
    axis([0 500 -100 100])
    
end



%% Plot Linear Filters switch
clear
date='20131025a'
path2load=['/Users/alexth/Desktop/old_stuff/ONresp/',date,'/'];
load([path2load,'LF_HC'])

path2save=[path2load,'/LF_plot_HC/'];
if ~exist(path2save,'dir')
    mkdir(path2save);
end



figure(1)
set(gcf,'Position', [1618 112 1224 825])
ndlist=unique(ndFullList,'stable');
for un=1:size(LinearFilter,3)
    
    for nd=1:length(ndlist)
        ind=find(ndFullList==ndlist(nd));
        a=ind(diff(ind)==1);
        
        subplot(2,3,nd)
        hold off
        plot(mean(LinearFilter(:,a(1):a(1)+1,un),2),'b','linewidth',2)
        hold on
        plot(mean(LinearFilter(:,a(2):a(2)+1,un),2),'r','linewidth',2)
        if length(a)==3
            plot(mean(LinearFilter(:,a(3):a(3)+1,un),2),'g','linewidth',2)
        end
        title(['ND',int2str(ndlist(nd))])
        
    end
    subplot(2,3,6)
    hold on
    plot(1,1,'b','linewidth',2)
    plot(1,1,'r','linewidth',2)
    plot(1,1,'g','linewidth',2)
    legend({'first time','after darker','after brighter'})
    
    subplot('Position',[0.5 0.95 0.00001 0.00001])
    set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
    title([names{un},'  HIGH contrast'],'FontSize',12,'FontWeight','bold','Interpreter','None')
    saveas(gcf,[path2save,names{un},'.png'])
    
end
