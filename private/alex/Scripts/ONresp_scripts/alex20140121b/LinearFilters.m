cd('/Users/alexth/Desktop/scripts/ONresp_scripts')

%% Prepare matrix for Full Field Flicker
clear
date='20140121b'
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

%correct!
ndFullList=[8 8 7 7 7 6 6 6 5 5 5 4 4 4];

save([path2save,'protocols_',codeWord],'ndFullList','mainpath','date','flicker','startTimes','maxLength','correctedProtocols','codeWord','file_list')



%% Get filter

clear
date='20140121b'
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


%% Plot Linear Filters
clear
date='20140121b'
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



% ND6 on 1 figure
figure(2)
set(gcf,'Position', [1 352 1011 625])
for un=1:size(LinearFilter,3)
    subplot(3,5,un)
    plot(LinearFilter(:,ndFullList==6,un))
    title(names{un},'interpreter','none')
    axis tight    
    line([0 500],[0 0],'color','k')    
end

