cd('/mnt/muench_data/user/alexandra/scripts')

%% Prepare matrix for Full Field Flicker
clear
cd('/mnt/muench_data/user/alexandra/scripts')
date='20121023'
codeWord='seq'

mainpath=['/mnt/muench_data/data/alexandra/MEA_data/',date,'/'];
hekapath=['/mnt/muench_data/data/alexandra/MEA_data/',date,'/HEKA/'];
heka=dir([hekapath,'*.phys']);

path2save=['/mnt/muench_data/data/alexandra/MEA_data/analysis/',date,'/'];
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

for i=1:length(file_list)
    protocol=read_header_field_heka(hekapath, heka(file_list(i)).name, 'Stimulus Protocol');
    protocol(protocol(:,3)==0&protocol(:,4)==0,4)=-1;
    protocol=round(protocol(2:end-2,[1,4]));
    protocol(:,1)=protocol(:,1)+tmp;
    correctedProtocols(:,:,i)=round(protocol);
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

save([path2save,'protocols_',codeWord],'mainpath','date','flicker','startTimes','maxLength','correctedProtocols','codeWord','file_list')



%% Get partial filter

clear
date='20121120_1'
codeWord='seq'
path2save=['/mnt/muench_data/data/alexandra/MEA_data/analysis/',date,'/'];
load([path2save,'protocols_',codeWord])

per=3600; % period in frames (start period). For sine:130, for HL10: 600, for H30s and L30s: 1800, for ffflicker: 1800, for 1 minute: 3600
dur=per*1; % duration of interval. for sine: per*5, for others: per*1

dim=floor((size(correctedProtocols,1)-1)/per);
while (size(correctedProtocols,1)-1)<(per*(dim-1)+dur)
    dim=dim-1;
end
filter_length=500;


units=dir([mainpath, 'units/*.mat']);

names=cell(1,length(units));
LinearFilter=zeros(filter_length,dim,length(file_list),length(units));
SpikeCount=zeros(length(file_list),length(units),dim);
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
            convSpikes=convolved(spikes,40,size(flicker,1)); 
            convSpikes=convSpikes(121:end-120);
            spikes=spikes-startTimes(i)+1;
            for t=1:dim
                startPoint=correctedProtocols(per*(t-1)+1,1,i)+1;
                endPoint=correctedProtocols(per*(t-1)+dur,1,i);
                spikesTMP=spikes(spikes>startPoint+500&spikes<endPoint);
                
                n=zeros(length(spikesTMP),filter_length);
                for k=1:filter_length
                    n(:,k)=flicker(spikesTMP-k+1,i);
                end
                LinearFilter(1:filter_length,t,i,cnt)=sum(n);
                SpikeCount(i,cnt,t)=length(spikesTMP);
                
            end
        end

        
        if isempty(regexp(units(cnt).name,date, 'once'))
            names{cnt}=[units(cnt).name(1),date,units(cnt).name(10:end-4)];
        else
            names{cnt}=units(cnt).name(1:end-4);
        end
    end
end
save([path2save,'LF_',codeWord],'LinearFilter','SpikeCount','names','units','per','dur','dim')


%% Quick
clear
date='20121120_1'
codeWord='quick';

mainpath=['/mnt/muench_data/data/alexandra/MEA_data/',date,'/'];
hekapath=['/mnt/muench_data/data/alexandra/MEA_data/',date,'/HEKA/'];
heka=dir([hekapath,'*.phys']);

path2save=['/mnt/muench_data/data/alexandra/MEA_data/analysis/',date,'/'];
if ~exist(path2save,'dir')
    mkdir(path2save);
end

% select heka files
file_list=[];
ff=[];
for i=1:length(heka)
    if ~isempty(regexp(heka(i).name,codeWord, 'once'))
        file_list=[file_list i];
    end
    if ~isempty(regexp(heka(i).name,'nd_8', 'once'))
        ff=[ff i];
    end
end
if ~isempty(ff)
    file_list(file_list<ff(1))=[];
end

% get protocol times
protocol=read_header_field_heka(hekapath, heka(file_list(1)).name, 'Stimulus Protocol');
protocol=round(protocol(2:end-2,[1,2]));
protocol(protocol(1:end,2)==1,2)=50;
protocol(protocol(1:end,2)==0,2)=10;
protocol(protocol(1:end,2)==-1,2)=30;

units=dir([mainpath,'units/*.mat']);
white_flash=zeros(4500,length(file_list),length(units));
black_flash=zeros(4500,length(file_list),length(units));
names=cell(length(units),1);
for cnt=1:length(units)
    
    load([mainpath,'units/',units(cnt).name]);
    
    for i=1:length(file_list)
        spikes=round(cell2mat(unit{1,2}(file_list(i),2))*1000/cell2mat(unit{1,1}(5,2))); %now in ms
        conv=convolved(spikes,40,75000);
        conv=conv(121:end-120);
        for st=1:4:20
            white_flash(:,i,cnt)=white_flash(:,i,cnt)+conv(protocol(st,1)-499:protocol(st,1)+4000)'/5;
            black_flash(:,i,cnt)=black_flash(:,i,cnt)+conv(protocol(st+2,1)-499:protocol(st+2,1)+4000)'/5;
        end
    end  
    
    if isempty(regexp(units(cnt).name,date, 'once'))
        names{cnt}=[units(cnt).name(1),date,units(cnt).name(10:end-4)];
    else
        names{cnt}=units(cnt).name(1:end-4);
    end
end

save([path2save,'quick'],'white_flash','black_flash','names')

%% Check polarity
date='20121120_1'
codeWord='seq'
path2save=['/mnt/muench_data/data/alexandra/MEA_data/analysis/',date,'/'];
load([path2save,'LF_',codeWord])
% find cells with stable filters at each ND
figure
for cnt=1:size(LinearFilter,4)
    subplot(4,6,cnt)
    hold on
    a=reshape(LinearFilter(:,1,1:24:floor(size(LinearFilter,3)/24)*24,cnt),500,floor(size(LinearFilter,3)/24));
    plot(a)
end 


%% zero crossing
% 12 trials
clear
boundaries=[155 130 115 85 75 75 75 75 85 85];
date='20121120_1';
codeWord='seq'
path2save=['/mnt/muench_data/data/alexandra/MEA_data/analysis/',date,'/'];
load([path2save,'protocols_',codeWord])
load([path2save,'LF_',codeWord])
% onOff=[1 1 1 -1 0 0 1 1 1 -1 1 -1 1 0 -1 1 -1 1 1 1 -1 1 1 1 -1 -1 1 1 -1 1 1 -1 0 1 1 1 0 1 1 1];%20121120
 onOff=[0 0 0 0 0 -1 -1 -1 1 0 1 1 1 1 -1 1 -1 0 0 1 -1 1 -1 -1]%20121120_1
cnt=1;
for i=1:24:24*7
    maxx=min(i+23,size(SpikeCount,1));
    for j=1:size(LinearFilter,4)
        t=reshape(LinearFilter(:,1,i:2:maxx,j),500,length(i:2:maxx));
        t=t-repmat(mean(t),size(t,1),1);
        t=t./repmat(sum(abs(t)),size(t,1),1);
        m=nanmean(t,2);
        if std(m(50:150))>std(m(350:450))*2
            if onOff(j)<0
                zc_HC(j,cnt)=find(m(boundaries(cnt):end)>0,1)+boundaries(cnt)-1;
            else
                zc_HC(j,cnt)=find(m(boundaries(cnt):end)<0,1)+boundaries(cnt)-1;
            end
            peak_HC(j,cnt)=myFit(m(zc_HC(j,cnt)-boundaries(cnt)+1:zc_HC(j,cnt)),zc_HC(j,cnt)-boundaries(cnt)+1,onOff(j),boundaries(cnt)-1);
        else
            peak_HC(j,cnt)=0;
            zc_HC(j,cnt)=0;
        end
        
        t=reshape(LinearFilter(:,1,i+1:2:maxx,j),500,length(i:2:maxx));
        t=t-repmat(mean(t),size(t,1),1);
        t=t./repmat(sum(abs(t)),size(t,1),1);
        m=nanmean(t,2);
        if std(m(50:150))>std(m(350:450))*2
            if onOff(j)<0
                zc_LC(j,cnt)=find(m(boundaries(cnt):end)>0,1)+boundaries(cnt)-1;
            else
                zc_LC(j,cnt)=find(m(boundaries(cnt):end)<0,1)+boundaries(cnt)-1;
            end
            peak_LC(j,cnt)=myFit(m(zc_LC(j,cnt)-boundaries(cnt)+1:zc_LC(j,cnt)),zc_LC(j,cnt)-boundaries(cnt)+1,onOff(j),boundaries(cnt)-1);
        else
            peak_LC(j,cnt)=0;
            zc_LC(j,cnt)=0;
        end
    end
    cnt=cnt+1;
end

save(['/mnt/muench_data/data/alexandra/MEA_data/analysis/zeroCrossing/zeroCrossing_',date],'zc_HC','zc_LC','peak_HC','peak_LC')

%%
clear
boundaries=[155 130 115 85 75 75 75 75 85 85];
date='20121004';
codeWord='seq'
path2save=['/mnt/muench_data/data/alexandra/MEA_data/analysis/',date,'/'];
load([path2save,'protocols_',codeWord])
load([path2save,'LF_',codeWord])
onOff=[1 -1 1 -1 1 -1 -1 -1 -1 1 -1 -1 -1 1 0 -1 -1 1 -1 -1 1 1 0 1 -1 -1 -1 -1];% 20121004

start=[27 36 45 54 63 72 73 97 121 145];
take=[9 9 9 9 9 1 24 24 24 24];
lc=[1 1 1 1 1 1 2 2 2 2];
cnt=1;
for ptr=1:length(take)
    i=start(ptr);

    for j=1:size(LinearFilter,4)
        t=reshape(LinearFilter(:,1,i:lc(ptr):i+take(ptr)-1,j),500,take(ptr)/lc(ptr));
        t=t-repmat(mean(t),size(t,1),1);
        t=t./repmat(sum(abs(t)),size(t,1),1);
        m=nanmean(t,2);
        if std(m(50:150))>std(m(350:450))*2
            if onOff(j)<0
                zc_HC(j,cnt)=find(m(boundaries(cnt):end)>0,1)+boundaries(cnt)-1;
            else
                zc_HC(j,cnt)=find(m(boundaries(cnt):end)<0,1)+boundaries(cnt)-1;
            end
            peak_HC(j,cnt)=myFit(m(zc_HC(j,cnt)-boundaries(cnt)+1:zc_HC(j,cnt)),zc_HC(j,cnt)-boundaries(cnt)+1,onOff(j),boundaries(cnt)-1);
        else
            peak_HC(j,cnt)=0;
            zc_HC(j,cnt)=0;
        end
        if lc(ptr)>1
            t=reshape(LinearFilter(:,1,(i+1):lc(ptr):i+take(ptr)-1,j),500,take(ptr)/lc(ptr));
            t=t-repmat(mean(t),size(t,1),1);
            t=t./repmat(sum(abs(t)),size(t,1),1);
            m=nanmean(t,2);
            if std(m(50:150))>std(m(350:450))*2
                if onOff(j)<0
                    zc_LC(j,cnt)=find(m(boundaries(cnt):end)>0,1)+boundaries(cnt)-1;
                else
                    zc_LC(j,cnt)=find(m(boundaries(cnt):end)<0,1)+boundaries(cnt)-1;
                end
                peak_LC(j,cnt)=myFit(m(zc_LC(j,cnt)-boundaries(cnt)+1:zc_LC(j,cnt)),zc_LC(j,cnt)-boundaries(cnt)+1,onOff(j),boundaries(cnt)-1);
            else
                peak_LC(j,cnt)=0;
                zc_LC(j,cnt)=0;
            end
        end
    end
    cnt=cnt+1;
end

save(['/mnt/muench_data/data/alexandra/MEA_data/analysis/zeroCrossing/zeroCrossing_',date],'zc_HC','zc_LC','peak_HC','peak_LC')

%% firing rate
% 12 trials
cd('/mnt/muench_data/user/alexandra/scripts')
clear
date='20121120_1'
codeWord='seq'
path2save=['/mnt/muench_data/data/alexandra/MEA_data/analysis/',date,'/'];
load([path2save,'protocols_',codeWord])
load([path2save,'LF_',codeWord])

% load(['/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_',date])

FiringRate=zeros(65000,size(SpikeCount,1),size(SpikeCount,2));
for cnt=1:length(units)
    cnt
    load([mainpath,'units/',units(cnt).name]);
    for i=1:length(file_list)
        spikes=round(cell2mat(unit{1,2}(file_list(i),2))*1000/cell2mat(unit{1,1}(5,2))); %now in ms
        if ~isempty(spikes)
            tmp=convolved(spikes,40,65000);
            FiringRate(:,i,cnt)=tmp(121:end-120);
        end
    end
end


for cnt=1:size(FiringRate,3)
    i=1;
    maxx=min(i+23,size(SpikeCount,1));
    for km=1:24:24*7
        maxx=min(km+23,size(SpikeCount,1));
        tmp=FiringRate(:,km:2:maxx,cnt);
        tmp=reshape(tmp(2501:62500,:),60000*length(km:2:maxx),1);
        
        frSTD_HC(i,cnt)=std(tmp);
        frMean_HC(i,cnt)=mean(tmp);
        
        tmp=FiringRate(:,(km+1):2:maxx,cnt);
        tmp=reshape(tmp(2501:62500,:),60000*length(km:2:maxx),1);
        frSTD_LC(i,cnt)=std(tmp);
        frMean_LC(i,cnt)=mean(tmp);
        a=FiringRate(200:1950,km:maxx,cnt);
        frspont(i,cnt)=mean(a(:));
        frSTDspont(i,cnt)=std(a(:));
        i=i+1;
    end
end

formula=1:7;
save(['/mnt/muench_data/data/alexandra/MEA_data/analysis/firingRateRatios/FiringRateRatio_',date],'frspont', 'frSTDspont','frSTD_HC','frMean_HC','frSTD_LC','frMean_LC','formula')


%% experiment summary
clear
date='20121120_1'
codeWord='seq'
path2save=['/mnt/muench_data/data/alexandra/MEA_data/analysis/',date,'/'];
load([path2save,'LF_',codeWord])
load([path2save,'protocols_',codeWord],'mainpath','file_list','correctedProtocols')
load(['/mnt/muench_data/data/alexandra/MEA_data/analysis/firingRateRatios/FiringRateRatio_',date],'frSTD_HC','frMean_HC','frSTD_LC','frMean_LC','frspont', 'frSTDspont','formula')
load(['/mnt/muench_data/data/alexandra/MEA_data/analysis/zeroCrossing/zeroCrossing_',date],'zc_HC','zc_LC','peak_HC','peak_LC')
% onOff=[-1 1 -1 -1 1 -1 0 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 0 1 0 0 -1 0 -1 0 1 1 -1 -1];% 20121002
% onOff=[-1 -1 1 -1 -1 1 1 1 -1 -1 1 -1 1 0 1 -1 1 -1 -1 -1 1 0 -1 1];% 20120928
% onOff=[1 -1 1 -1 1 -1 -1 -1 -1 1 -1 -1 -1 1 0 -1 -1 1 -1 -1 1 1 0 1 -1 -1 -1 -1];% 20121004
% onOff=[1 1 1 -1 0 0 1 1 1 -1 1 -1 1 0 -1 1 -1 1 1 1 -1 1 1 1 -1 -1 1 1 -1 1 1 -1 0 1 1 1 0 1 1 1];%20121120
 onOff=[0 0 0 0 0 -1 -1 -1 1 0 1 1 1 1 -1 1 -1 0 0 1 -1 1 -1 -1]%20121120_1
a=who;
str='';
for i=1:length(a)
    str=[str,', ''',a{i},''''];
end
eval(['save(','''/mnt/muench_data/data/alexandra/MEA_data/analysis/summary_',date,'''', str,')'])
