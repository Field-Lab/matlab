clear
addpath(genpath('/Users/alexth/Desktop/Scripts'))



%% parameters for filter

pathName='/Users/alexth/Desktop/myCNGa3/';
date='20120928';
contrastKey='seq';
params.per=3600; % frames to take for calculations
params.filter_length=500;


%% Prepare matrix for Full Field Flicker
unitspath=fullfile(pathName,date,'units',filesep);
hekapath=fullfile(pathName,date,'HEKA',filesep);
heka=dir([hekapath,'*.phys']);

path2save=fullfile(pathName,date,'analysis',filesep);
if ~exist(path2save,'dir')
    mkdir(path2save);
end

% select heka files
file_list=[];
for i=1:length(heka)
    if ~isempty(regexp(heka(i).name,contrastKey, 'once'))
        file_list=[file_list i];
    end
end

% get protocol correction matrix
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

save([path2save,'protocols_',contrastKey],'ndFullList','unitspath','date','flicker','startTimes','maxLength','correctedProtocols','contrastKey','file_list')



%% Get filter

units=dir([unitspath, '*.mat']);

names=cell(1,length(units));
LinearFilter=zeros(params.filter_length,length(file_list),length(units));
SpikeCount=zeros(length(file_list),length(units));
if length(unique(startTimes))>1
    for i=1:size(correctedProtocols,3)
        correctedProtocols(:,1,i)=correctedProtocols(:,1,i)-startTimes(i);
    end
else
    correctedProtocols(:,1,:)=correctedProtocols(:,1,:)-startTimes(1);
end

for cnt=1:length(units)
    
    load([unitspath,units(cnt).name]);
    for i=1:length(file_list)
        spikes=round(cell2mat(unit{1,2}(file_list(i),2))*1000/cell2mat(unit{1,1}(5,2))); %now in ms
        if ~isempty(spikes)
            spikes=spikes-startTimes(i)+1;
            startPoint=correctedProtocols(1,1,i)+1;
            endPoint=correctedProtocols(params.per,1,i);
            spikesTMP=spikes(spikes>startPoint+500&spikes<endPoint);
            
            n=zeros(length(spikesTMP),params.filter_length);
            for k=1:params.filter_length
                n(:,k)=flicker(spikesTMP-k+1,i);
            end
            LinearFilter(1:params.filter_length,i,cnt)=sum(n);
            SpikeCount(i,cnt)=length(spikesTMP);
        end

    end    
   
    if isempty(regexp(units(cnt).name,date, 'once'))
        names{cnt}=[units(cnt).name(1),date,units(cnt).name(10:end-4)];
    else
        names{cnt}=units(cnt).name(1:end-4);
    end
end

%polarity check
pc1=zeros(size(LinearFilter,3),size(LinearFilter,2));
for cnt=1:144
    tmp=squeeze(LinearFilter(1:200,cnt,:));
    x=tmp';
    x=cov(x);
    [V,~]=eig(x);
    pc_vectors=V(:,end);
    pc1(:,cnt)=tmp'*pc_vectors;
end
tmp=sum(pc1>0,2)>sum(pc1<0,2);
polarity=zeros(1,length(units));
polarity(tmp)=-1; %OFF cell
polarity(~tmp)=1; %ON cell

filtersFile=[path2save,'filters_',contrastKey];

save(filtersFile,'ndFullList','date','LinearFilter','SpikeCount','polarity','names','units','params','filtersFile')


%% fit
allthere=whos;
for i=1:size(allthere,1)
    if ~strcmp(allthere(i).name,'filtersFile')&&~strcmp(allthere(i).name,'allthere')
        eval(['clear ', allthere(i).name])
    end
end
clear allthere

global fitParams currentUnit start_trial end_trial currentLFs trialsPerScreen fitResults
load(filtersFile)

currentUnit=1;

savePath=filtersFile;
trialsPerScreen=min(size(LinearFilter,2),49);
pol=(polarity(currentUnit)/2)+0.5;
start_trial=1;
end_trial=trialsPerScreen;
nUnits=size(LinearFilter,3);
nTrials=size(LinearFilter,2);
unitName=names{currentUnit};
start_trial=1;
end_trial=min(nTrials,49);
currentLFs=LinearFilter(:,:,currentUnit);

fitResults=zeros(3,nTrials,nUnits);

handles=creatFigure(unitName,pol,currentUnit,nUnits,nTrials,start_trial,end_trial,trialsPerScreen);

calcParams(currentLFs,SpikeCount(:,currentUnit),handles,trialsPerScreen,start_trial);



