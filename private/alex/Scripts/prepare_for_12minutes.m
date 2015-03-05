cd('/mnt/muench_data/user/alexandra/scripts')

%% Prepare matrix for Full Field Flicker
clear
cd('/mnt/muench_data/user/alexandra/scripts')
date='20121023';
codeWord='seq'

mainpath=['/mnt/muench_data/data/alexandra/MEA_data/',date,'/'];
hekapath=['/mnt/muench_data/data/alexandra/MEA_data/',date,'/HEKA/'];
heka=dir([hekapath,'*.phys']);

if ~exist([mainpath 'LinearFilters/'],'dir')
    mkdir([mainpath, 'LinearFilters']);
end
path2save=[mainpath 'LinearFilters/'];

% make list of files to put into easy formatted unit (choosing by the name of the stimulus in heka files)
file_list=[];
for i=1:length(heka)
    if ~isempty(regexp(heka(i).name,codeWord, 'once'))
        file_list=[file_list i];
    end
end


%get list of units to be reformatted
basic_format_units_list=dir([mainpath, 'units/*.mat']);

% get protocol
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

save([path2save,date,'_',codeWord,'_protocols_raw'],'protocol_accum')

% load([path2save,date,'_',codeWord,'_protocols_raw'])

tmp=(0:size(protocol_accum,1))'*(7.85/3600);
tmp=repmat(tmp(1:end-1),1,size(protocol_accum,3));
abused_accum=reshape(protocol_accum(:,1,:),size(protocol_accum,1),length(file_list));
abused_accum=round(abused_accum+tmp);

flicker=zeros(abused_accum(end,1),length(file_list));
startingPoints=zeros(1,length(file_list));
maxLength=zeros(1,length(file_list));

for i=1:length(file_list)
    flips=abused_accum(1:end-2,i);
    pxls=(protocol_accum(1:end-2,2,i)-30)/30;
    startingPoints(i)=flips(1);
    flips=flips-flips(1)+1;
    
    tmp=zeros(flips(end),1);
    pxls(pxls==0)=100;
    tmp(flips(2:end),1)=pxls(2:end);
    % fill stimulus down to ms with brightness values where spikes
    % happen
    droppedMax=max(diff(flips));
    for j=1:droppedMax+1
        subFrame=flips(2:end)-j;
        subFrame(subFrame<1)=[];
        subFrame(find(tmp(subFrame,1)))=[];
        tmp(subFrame,1)=tmp(subFrame+1,1);
    end
    tmp(tmp(:,1)==100,1)=0;
    flicker(1:length(tmp),i)=tmp;
    maxLength(i)=length(tmp);
end
flicker=flicker(1:max(maxLength),:);

save([path2save,date,'_flicker_',codeWord],'flicker','startingPoints','maxLength','abused_accum')



%% Get partial filter

clear
dates=cell(1,3);
dates{1}='20121023';
dates{2}='20130301_1';
dates{3}='20130301_2';

codeWord='seq';
per=3600; % period in frames (start period). For sine:130, for HL10: 600, for H30s and L30s: 1800
dur=per*1; % duration of interval. for sine: per*5, for others: per*1

for ttt=dates(1)
    date=cell2mat(ttt)
    mainpath=['/mnt/muench_data/data/alexandra/MEA_data/',date,'/'];
    basic_format_units_list=dir([mainpath, 'units/*.mat']);
    path2take=[mainpath,'LinearFilters/'];
    load([path2take,date,'_flicker_',codeWord])
    dim=floor((size(abused_accum,1)-2)/per);
    while (size(abused_accum,1)-2)<(per*(dim-1)+dur)
        dim=dim-1;
    end
    filter_length=500;
    
    hekapath=['/mnt/muench_data/data/alexandra/MEA_data/',date,'/HEKA/'];
    heka=dir([hekapath,'*.phys']);
    file_list=[];
    for i=1:length(heka)
        if ~isempty(regexp(heka(i).name,codeWord, 'once'))
            file_list=[file_list i];
        end
    end
    
    path2save=['/mnt/muench_data/data/alexandra/MEA_data/',date,'/',codeWord,'_per_',int2str(per),'_dur_',int2str(dur),'/'];
    
    if ~exist(path2save,'dir')
        mkdir(path2save);
    end
    
    names=cell(1,length(basic_format_units_list));
    LinearFilter=zeros(filter_length,dim,length(file_list),length(basic_format_units_list));
    SpikeCount=zeros(length(file_list),length(basic_format_units_list),dim);
    
    for cnt=1:length(basic_format_units_list)
        
        load([mainpath,'units/',basic_format_units_list(cnt).name]);
        for i=1:length(file_list)
            spikes=round(cell2mat(unit{1,2}(file_list(i),2))*1000/cell2mat(unit{1,1}(5,2))); %now in ms
            if ~isempty(spikes)
                spikes=spikes-startingPoints(i)+1;
                spikes(spikes<filter_length|spikes>size(flicker,1))=[];
                for t=1:dim
                    spikesTMP=spikes(spikes>(abused_accum(per*(t-1)+1,i)-startingPoints(i))&spikes<(abused_accum(per*(t-1)+dur,i)-startingPoints(i)));
                    n=zeros(length(spikesTMP),filter_length);
                    for k=1:filter_length
                        n(:,k)=flicker(spikesTMP-k+1,i);
                    end
                    LinearFilter(1:filter_length,t,i,cnt)=sum(n);
                    SpikeCount(i,cnt,t)=length(spikesTMP);
                end
            end
            names{cnt}=basic_format_units_list(cnt).name;
        end
    end
    save([path2save,date,'_LF_',codeWord,'_',int2str(per),'_',int2str(dur)],'LinearFilter','SpikeCount','names')
    
end

%% Nonlinearity
clear
dates=cell(1,3);
dates{1}='20121023';
dates{2}='20130301_1';
dates{3}='20130301_2';

codeWord='seq';
per=3600; % period in frames (start period)
dur=per*1;

date='20121023';
path2load=['/mnt/muench_data/data/alexandra/MEA_data/',date,'/',codeWord,'_per_',int2str(per),'_dur_',int2str(dur),'/'];
load([path2load,date,'_LF_',codeWord,'_',int2str(per),'_',int2str(dur)])


for ttt=dates(1)
    date=cell2mat(ttt)
    mainpath=['/mnt/muench_data/data/alexandra/MEA_data/',date,'/'];
    basic_format_units_list=dir([mainpath, 'units/*.mat']);
    path2take=[mainpath,'/LinearFilters/'];
    load([path2take,date,'_flicker_',codeWord])
    dim=floor((size(abused_accum,1)-2)/per);
    while (size(abused_accum,1)-2)<(per*(dim-1)+dur)
        dim=dim-1;
    end
    
    path2load=['/mnt/muench_data/data/alexandra/MEA_data/',date,'/',codeWord,'_per_',int2str(per),'_dur_',int2str(dur),'/'];
    load([path2load,date,'_LF_',codeWord,'_',int2str(per),'_',int2str(dur)])
    
    
    hekapath=['/mnt/muench_data/data/alexandra/MEA_data/',date,'/HEKA/'];
    heka=dir([hekapath,'*.phys']);
    file_list=[];
    for i=1:length(heka)
        if ~isempty(regexp(heka(i).name,codeWord, 'once'))
            file_list=[file_list i];
        end
    end
    
    path2save=['/mnt/muench_data/data/alexandra/MEA_data/',date,'/',codeWord,'_per_',int2str(per),'_dur_',int2str(dur),'/'];
    
    if ~exist(path2save,'dir')
        mkdir(path2save);
    end
    
    nonLinear=zeros(4,dim,size(LinearFilter,3),length(basic_format_units_list));
    nonlin=zeros(51,1);
    
    convSpikes=zeros(size(flicker,1),size(LinearFilter,3),length(basic_format_units_list));
    for cnt=1:length(basic_format_units_list)
        load([mainpath,'units/',basic_format_units_list(cnt).name]);
        for i=1:length(file_list)
            spikes=round(cell2mat(unit{1,2}(file_list(i),2))*1000/cell2mat(unit{1,1}(5,2))); %now in ms
            spikes=spikes-startingPoints(i)+1;
            spikes(spikes<500|spikes>size(flicker,1))=[];
            tmp=convolved(spikes,40,size(flicker,1));
            convSpikes(:,i,cnt)=tmp(121:end-120);
        end
    end
    
    for i=1:length(file_list)
        i
        tic
        for t=1:dim
            subflicker=flicker(abused_accum(per*(t-1)+1,i)-startingPoints(i)+1:abused_accum(per*(t-1)+dur,i)-startingPoints(i),i);
            tmp=zeros(length(subflicker)-499,length(basic_format_units_list));
            a=reshape(LinearFilter(500:-1:1,t,i,:),500,length(basic_format_units_list))';
            
            for rt=500:length(subflicker)
                tmp(rt-499,:)=a*subflicker(rt-499:rt);
            end
            
            a=convSpikes((abused_accum(per*(t-1)+1,i)-startingPoints(i))+500:(abused_accum(per*(t-1)+dur,i)-startingPoints(i)),i,:);
            a=reshape(a,size(a,1),size(a,3));
            
            for cnt=1:length(basic_format_units_list)
                if sum(abs(tmp(:,cnt)))>5
                    stepSize=(max(tmp(:,cnt))-min(tmp(:,cnt)))/50;
                    x=min(tmp(:,cnt)):stepSize:max(tmp(:,cnt));
                    kk=1;
                    for rt=x
                        z=find(tmp(:,cnt)>=rt&tmp(:,cnt)<rt+stepSize);
                        nonlin(kk)=mean(a(z,cnt)');
                        kk=kk+1;
                    end
                    
                    while sum(isnan(nonlin))
                        nonlin(isnan(nonlin))=nonlin(isnan(nonlin(2:end)));
                    end
                    try
                        res=fit(x',nonlin,'exp1');
                        nonLinear(1,t,i,cnt)=res.a;
                        nonLinear(2,t,i,cnt)=res.b;
                        nonLinear(3,t,i,cnt)=x(1);
                        nonLinear(4,t,i,cnt)=x(end);
                    catch
                        disp('could NOT fit')
                        date
                        cnt
                        t
                        disp('END OF could NOT fit')
                    end
                    
                else
                    disp('TOO LOW FR')
                    i
                    t
                    cnt
                end
            end
        end
        toc
    end
    
    save([path2save,date,'_LF_',codeWord,'_',int2str(per),'_',int2str(dur)],'LinearFilter','nonLinear','SpikeCount','names')
    
    
end


%% Accumulate data
clear
codeWord='seq';
per=3600;
dur=per*1;

date='20121023';
path2load=['/mnt/muench_data/data/alexandra/MEA_data/',date,'/',codeWord,'_per_',int2str(per),'_dur_',int2str(dur),'/'];
load([path2load,date,'_LF_',codeWord,'_',int2str(per),'_',int2str(dur)])

% sigALL: contrast function (272,6 or 1 points)
path2load=['/mnt/muench_data/data/alexandra/MEA_data/',date,'/LinearFilters/'];
load([path2load,date,'_',codeWord,'_protocols_raw'])
sigALL=zeros(1,size(LinearFilter,3));
for k=1:size(LinearFilter,3) 
    protocol=protocol_accum(:,2,k);
    sigALL(k)=std(protocol((i-1)*per+1:(i-1)*per+dur));
end
% % average contrast (58 points)
% t=reshape(sigALL(41-13:end-13),periodsTake,floor(272/periodsTake));
% contrChange=mean(t');

clear abused_accum protocol_accum protocol 

% cell polarity by filter at ND6
onOff=zeros(1,size(LinearFilter,4))-1;
col='';
minimax=onOff;
for cellID=1:size(LinearFilter,4)    
    tmp=LinearFilter(80:180,:,27,cellID);
    [~, b]=min(tmp);
    [~, c]=max(tmp);
    if b>c % on cell
        onOff(cellID)=1;
        minimax(cellID)=c+79;
        col(cellID)='r';
    else
        minimax(cellID)=b+79;
        col(cellID)='b';
    end
end


% index of zero crossing
clear ind
for cellID=1:size(LinearFilter,4)    
    for k=1:size(LinearFilter,3) 
        tmp=LinearFilter(:,:,k,cellID);        
        for j=1:dim
            b=tmp(:,j);
            if onOff(cellID)>0
                [~, minPos]=max(b(1:150));
            else
                [~, minPos]=min(b(1:150));
            end
            temp=find(onOff(cellID)*b(minPos:end)<0,1)+minPos-1;
            if isempty(temp)
                ind(j,k,cellID)=0;
            else
                ind(j,k,cellID)=find(onOff(cellID)*b(minPos:end)<0,1)+minPos-1;
            end
        end
    end
end

% peak value of the filter
peak=zeros(dim,size(LinearFilter,3),size(LinearFilter,4));
peakLatency=zeros(dim,size(LinearFilter,3),size(LinearFilter,4));
for k=1:size(LinearFilter,3) 
    for cellID=1:size(LinearFilter,4)  
        if sum(ind(:,k,cellID))>0
            t=LinearFilter(:,:,k,cellID);
            if onOff(cellID)>0
                [peak(:,k,cellID),peakLatency(:,k,cellID)]=max(t(1:max(ind(:,k,cellID)),:));
            else
                [peak(:,k,cellID),peakLatency(:,k,cellID)]=min(t(1:max(ind(:,k,cellID)),:));
            end
        end
    end
end

% adjusted gain
nonLinear_adjusted=zeros(dim,size(LinearFilter,3),size(LinearFilter,4));
for k=1:size(LinearFilter,3)
    for cellID=1:size(LinearFilter,4)
        if sum(ind(:,k,cellID))>0
            t=LinearFilter(:,:,k,cellID);
            if onOff(cellID)>0
                p=max(t(1:max(ind(:,k,cellID)),:));
            else
                p=min(t(1:max(ind(:,k,cellID)),:));
            end
            nonLinear_adjusted(:,k,cellID)=nonLinear(2,:,k,cellID)./abs(p);
        end
    end
end


save(['/mnt/muench_data/data/alexandra/MEA_data/contrast12min_analysis/',codeWord,'_accumParameters'],'names','LinearFilter','SpikeCount',...
    'nonLinear','sigALL','onOff','ind','peak','peakLatency','nonLinear_adjusted','col','dim','per','dur')


%% Plots
cd('/mnt/muench_data/user/alexandra/scripts')

clear
codeWord='seq';
load(['/mnt/muench_data/data/alexandra/MEA_data/contrast12min_analysis/',codeWord,'_accumParameters'])


% firing rate at high vs firing rate at low
figure
nds='876543212345';
ranges=[0 60 0 60];
cnt=1;
for i=1:24:24*12
    subplot(4,3,cnt)
    hold on
    c=onOff<0;
    a=mean(SpikeCount(i:2:i+23,c))/60; % roughly Hz
    b=mean(SpikeCount(i+1:2:i+23,c))/60;
    plot(a,b,'.b')
    [p,hoff(i)]=ttest(a-b);
    
    c=onOff>0;
    a=mean(SpikeCount(i:2:i+23,c))/60; % roughly Hz
    b=mean(SpikeCount(i+1:2:i+23,c))/60;
    plot(a,b,'.r')
    [p,hon(i)]=ttest(a-b);    
    
    xlabel('spike count at high')
    ylabel('spike count at low')

    axis(ranges)
    line(ranges(1:2),[0 0],'color','k')
    line([0 0],ranges(3:4),'color','k')
    title(['ND', nds(cnt),', red ON, blue OFF'])
    line(ranges(1:2),ranges(3:4),'color','k')
    line(ranges(1:2)+3,ranges(3:4),'color','k')
    line(ranges(1:2)-3,ranges(3:4),'color','k')
    cnt=cnt+1;
end 


% zero crossing at high vs zero crossing  at low
figure
nds='876543212345';
ranges=[0 300 0 300];
cnt=1;
for i=1:24:24*12
    subplot(4,3,cnt)
    hold on
    c=onOff<0;
    a=mean(ind(1,i:2:i+23,c));
    b=mean(ind(1,i+1:2:i+23,c));
    a=reshape(a,size(a,3),1);
    b=reshape(b,size(b,3),1);
    plot(a,b,'.b')
    [p,hoff(cnt)]=ttest(a-b);
    
    c=onOff>0;
    a=mean(ind(1,i:2:i+23,c));
    b=mean(ind(1,i+1:2:i+23,c));
    a=reshape(a,size(a,3),1);
    b=reshape(b,size(b,3),1);
    plot(a,b,'.r')
    [p,hon(cnt)]=ttest(a-b);    
    
    xlabel('zero x at high')
    ylabel('zero x at low')

    axis(ranges)
    line(ranges(1:2),[0 0],'color','k')
    line([0 0],ranges(3:4),'color','k')
    title(['ND', nds(cnt),', red ON, blue OFF'])
    line(ranges(1:2),ranges(3:4),'color','k')
    cnt=cnt+1;
end 


%firing rate vs zero crossing
figure
nds='876543212345';
ranges=[-10 10 -50 50]*2;
cnt=1;
for i=1:24:24*12
    subplot(4,3,cnt)
    hold on
    c=onOff<0;
    a=mean(SpikeCount(i:2:i+23,c))/60-mean(SpikeCount(i+1:2:i+23,c))/60;  
    b=mean(ind(1,i:2:i+23,c))-mean(ind(1,i+1:2:i+23,c));
    b=reshape(b,1,size(b,3));
    plot(a,b,'.b')    

    c=onOff>0;
    a=mean(SpikeCount(i:2:i+23,c))/60-mean(SpikeCount(i+1:2:i+23,c))/60;  
    b=mean(ind(1,i:2:i+23,c))-mean(ind(1,i+1:2:i+23,c));
    b=reshape(b,1,size(b,3));
    plot(a,b,'.r')
    
    
    xlabel('spike count high-low')
    ylabel('zero crossing high-low')
    axis(ranges)
    line(ranges(1:2),[0 0],'color','k')
    line([0 0],ranges(3:4),'color','k')
    title(['ND', nds(cnt),', red ON, blue OFF'])
    cnt=cnt+1;
end 

