cd('/Users/alexth/Desktop/scripts/stefano_scripts')

%% Prepare matrix for Full Field Flicker
clear
date='20130820'
codeWord='HC'

mainpath=['/Users/alexth/Desktop/old_stuff/',date,'/'];
hekapath=['/Users/alexth/Desktop/old_stuff/',date,'/HEKA/'];
heka=dir([hekapath,'*.phys']);

path2save=['/Users/alexth/Desktop/old_stuff/stefano_analysis/',date,'/'];
if ~exist(path2save,'dir')
    mkdir(path2save);
end

% % select heka files
% file_list=[];
% for i=1:length(heka)
%     if ~isempty(regexp(heka(i).name,codeWord, 'once'))
%         file_list=[file_list i];
%     end
% end

% select heka files from units
units=dir([mainpath, 'units/*.mat']);
load([mainpath,'units/',units(1).name]);
file_list=[];
cnt=1;
for i=1:size(unit{1,2},1)
    if ~isempty(regexp(unit{1,2}{i,1},codeWord, 'once'))
        heka_list{cnt}=unit{1,2}{i,1};
        file_list=[file_list i];
        cnt=cnt+1;
    end
end


% get protocol
protocol=read_header_field_heka(hekapath, heka_list{2}, 'Stimulus Protocol');
protocol(protocol(:,3)==0&protocol(:,4)==0,4)=-1;
protocol=round(protocol(2:end-2,[1,4]));
correctedProtocols=zeros(size(protocol,1),2,size(heka_list,2));
tmp=(0:size(correctedProtocols,1)-1)'*(7.85/3600);
ndFullList=zeros(1,size(heka_list,2));

for i=1:size(heka_list,2)
    hekaFile=heka_list{i};
    protocol=read_header_field_heka(hekapath, hekaFile, 'Stimulus Protocol');
    protocol(protocol(:,3)==0&protocol(:,4)==0,4)=-1;
    protocol=round(protocol(2:end-2,[1,4]));
    protocol(:,1)=protocol(:,1)+tmp;
    correctedProtocols(:,:,i)=round(protocol);
    
    a=regexp(hekaFile,'ND');
    a=hekaFile(a+2);
    if length(a)>1
        ndFullList(i)=str2num(a(1))+str2num(a(2));
    else
        ndFullList(i)=str2num(a(1));
    end
end
clear protocol

startTimes=reshape(correctedProtocols(1,1,:),1,size(heka_list,2));
endTimes=reshape(correctedProtocols(end,1,:),1,size(heka_list,2));
maxLength=max(endTimes-startTimes);
flicker=zeros(maxLength,size(heka_list,2));

for i=1:size(heka_list,2)
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


save([path2save,'protocols_',codeWord],'ndFullList','file_list','mainpath','date','flicker','startTimes','maxLength','correctedProtocols','codeWord','heka_list')



%% Get filter

clear
date='20130820'
codeWord='HC'

path2save=['/Users/alexth/Desktop/old_stuff/stefano_analysis/',date,'/'];
load([path2save,'protocols_',codeWord])

per=1800; % frames to take for calculations
filter_length=500;

units=dir([mainpath, 'units/*.mat']);

names=cell(1,length(units));
LinearFilter=zeros(filter_length,size(heka_list,2),length(units));
SpikeCount=zeros(size(heka_list,2),length(units));
if length(unique(startTimes))>1
    for i=1:size(correctedProtocols,3)
        correctedProtocols(:,1,i)=correctedProtocols(:,1,i)-startTimes(i);
    end
else
    correctedProtocols(:,1,:)=correctedProtocols(:,1,:)-startTimes(1);
end

for cnt=1:length(units)
    
    load([mainpath,'units/',units(cnt).name]);
    for i=1:size(heka_list,2)
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

save([path2save,'LF_',codeWord],'ndFullList','heka_list','file_list','LinearFilter','SpikeCount','names','units','per')


%% Plot Linear Filters
clear
date='20130820'

path2load=['/Users/alexth/Desktop/old_stuff/stefano_analysis/',date,'/'];
load([path2load,'LF_HC'])

path2save=[path2load,'/LF_plot_HC/'];
if ~exist(path2save,'dir')
    mkdir(path2save);
end


figure(1)
set(gcf,'Position', [1 352 1011 625])
nds='876';
for un=1:size(LinearFilter,3)
    for cnt=1:3
        subplot(2,2,cnt)
        a=LinearFilter(:,ndFullList==str2num(nds(cnt)),un);
        if size(a,2)>20
            hold off
            plot(a(:,1:end-10),'b')
            hold on
            plot(a(:,end-9:end),'r')
        else
            hold off
            plot(a)
        end
        title(['ND',nds(cnt)])
    end
    for i=1:3
        subplot(2,2,i)
        line([0 500],[0 0],'color','k')
        k=LinearFilter(:,:,un);
        axis([0 500 min(k(:)) max(k(:))])
    end
    subplot('Position',[0.5 0.95 0.00001 0.00001])
    set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
    title([names{un},'  HIGH contrast'],'FontSize',12,'FontWeight','bold','Interpreter','None')
    saveas(gcf,[path2save,names{un},'.png'])    
end


date='20130729'
figure
path2load=['/Users/alexth/Desktop/old_stuff/stefano_analysis/',date,'/'];
load([path2load,'LF_HC'])
for i=1:size(SpikeCount,2)
    subplot(7,8,i)
    plot(SpikeCount(:,i))
    hold on
end
path2load=['/Users/alexth/Desktop/old_stuff/stefano_analysis/',date,'/'];
load([path2load,'LF_LC'])



date='20130820'
codeWord='HC'
path2save=['/Users/alexth/Desktop/old_stuff/stefano_analysis/',date,'/'];
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
% 20130628b
% onOff=[1 -1 1 1 -1 0 0 0 -1 -1 -1 1 1 1 -1 1 1 -1 1 -1 1 -1 1 0 0 -1 -1 0 1 1 1 -1 0 0 -1 -1 1 -1 1 -1 1 -1 -1 -1 1 -1 0 -1 1 -1 1 -1];

% 20130629
onOff=[1 -1 1 -1 0 1 1 0 1 -1 1 1 1 -1 -1 1 0 0 0 1 1 -1 0 1 0 0 -1 0 1 0 0 0 0 0 0 0 0 0 -1 0 0 0 -1 -1 0 0 -1 1 -1 -1];


% 20130729
onOff=[-1 0 0 0 1 1 0 0 -1 0 -1 -1 -1 0 0 0 1 0 -1 1 -1 0 0 -1 -1 1 0 0 -1 0 0 0 0];

% 20130808
onOff=[1 0 0 1 -1 -1 1 -1 -1 1 1 -1 0 0 1 0 0 0 1 0 -1 1 1 -1 0 0 0 -1 0 0 1 0 0 0 1 -1 0 0 0 -1 1 1 0 -1 -1 1 1 0 1];


nds='876';
for cnt=1:size(FiringRate,3)
    for i=1:3
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
figure
nds='876';
for i=1:3
    subplot(1,3,i)
    plot(frspont(i,onOff>0),frMean_HC(i,onOff>0)-frspont(i,onOff>0),'*r')
    hold on
    plot(frspont(i,onOff<0),frMean_HC(i,onOff<0)-frspont(i,onOff<0),'*b')
    title(['ND',nds(i)])
    line([0 45],[0 0],'color','k')
end






nds='8765';
for cnt=1:size(FiringRate,3)
    for i=1:4
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



figure
nds='8765';
for i=1:4
    subplot(2,2,i)
    plot(frspont(i,onOff>0),frMean_HC(i,onOff>0)-frspont(i,onOff>0),'*r')
    hold on
    plot(frspont(i,onOff<0),frMean_HC(i,onOff<0)-frspont(i,onOff<0),'*b')
    title(['ND',nds(i)])
    line([0 45],[0 0],'color','k')
end


figure
nds='7655'
for i=1:4
    subplot(2,2,i)
    plot(frMean_HC(i,onOff>0)-frspont(i,onOff>0),'*r')
    hold on
    plot(frMean_HC(i,onOff<0)-frspont(i,onOff<0),'*b')
    title(['ND',nds(i)])
    line([0 22],[0 0],'color','k')
end

figure
for i=1:size(SpikeCount,2)
    subplot(7,8,i)
    plot(SpikeCount(:,i),'r')
    title(names{i},'Interpreter','None')

end



date='20130628b'
figure
path2load=['/gpfs01/muench/data/alexandra/MEA_data/stefano_analysis/',date,'/'];
load([path2load,'LF_HC'])
for i=1:size(SpikeCount,2)
    subplot(7,8,i)
    plot(mean(LinearFilter(:,24:73,i),2))
    hold on
     plot(mean(LinearFilter(:,74:123,i),2),'g')
     plot(mean(LinearFilter(:,124:173,i),2),'r')
     plot(mean(LinearFilter(:,174:223,i),2),'m')
end






figure(1)
set(gcf,'Position', [1 352 1011 625])
nds='876';
for un=1:33
%     for cnt=1:3
%         subplot(2,2,cnt)
%         plot(LinearFilter(:,ndFullList==str2num(nds(cnt)),un))
%         title(['ND',nds(cnt)])
%     end
%     

    k=50;
    for i=1:5
        tmp(1:500,i)=mean(LinearFilter(:,(24:27)+(i-1)*10+k,un),2);
    end
    subplot(6,6,un)
    plot(tmp(:,[1:5]))
    axis([0 500 -100 100])
%         
%     for i=1:3
%         subplot(2,2,i)
%         line([0 500],[0 0],'color','k')
%         k=LinearFilter(:,:,un);
%         axis([0 500 min(k(:)) max(k(:))])
%     end
%     subplot('Position',[0.5 0.95 0.00001 0.00001])
%     set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
%     title([names{un},'  HIGH contrast'],'FontSize',12,'FontWeight','bold','Interpreter','None')
%     saveas(gcf,[path2save,names{un},'.png'])
    
end
