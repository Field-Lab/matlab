clear

cd('/mnt/muench_data/user/alexandra/scripts')
date='20121207';

length2take=60; % time in s

path2data=['/mnt/muench_data/data/alexandra/MEA_data/',date,'/LFP/'];
keyWord='HCseq';
lfpData=dir([path2data,'*',keyWord,'*.mat']);

hekapath=['/mnt/muench_data/data/alexandra/MEA_data/',date,'/HEKA/'];
heka=dir([hekapath,'*',keyWord,'*.phys']);

path2save=['/mnt/muench_data/data/alexandra/MEA_data/',date,'/LFP_filters/'];
if exist(path2save,'dir')    
    rmdir(path2save,'s')
end
mkdir(path2save);


filters=zeros(500,length(lfpData),60);
k=zeros(length2take*1000-1000,500);
tic
for i=1:length(lfpData)
    i
    load([path2data,lfpData(i).name])
    protocol=read_header_field_heka(hekapath, heka(i).name, 'Stimulus Protocol');
    protocol(protocol(:,3)==0&protocol(:,4)==0,4)=-1;
    protocol=round(protocol(2:end,[1,4]));
    
    
    
    flips=protocol(:,1); % flips in ms, rounded
    flips(flips>length2take*1000)=[];
    pxls=((protocol(:,2)-30)/30); % pxs in ms, normalized -1:1
    pxls(length(flips)+1:end)=[];
    
    %         a=lfp(:,trial_counter);
    
    tmp=zeros(flips(end)+1,1); % 1 column - color, 2 column - spikes
    tmp([1; flips(1:end-1)])=pxls;
    
    % fill stimulus down to ms with brightness values where spikes
    % happen
    for j=1:50
        subscr=flips(1:end-1)+j;
        while length(tmp)<subscr(end)
            subscr(end)=[];
        end
        subscr(tmp(subscr,1)~=0)=[];
        tmp(subscr,1)=tmp(subscr-1,1);
    end
    
    for j=1:length2take*1000-1000
        k(j,:)=tmp(j+499:-1:j);
    end
    
    for j=1:60
        filters(:,i,j)=lfp(501:length2take*1000-1000+500,j)'*k;
    end
    
end
toc
save([path2save,'lfp_filters_',keyWord],'filters')




figure
for i=1:60
    subplot(8,8,i)
    plot(filters(:,30,i)/5e6)
    title(int2str(i))
    axis([1 500 -1 0.5])
end
subplot(8,8,61)
plot(mean(filters(:,30,:)/5e6,3))
axis([1 500 -1 0.5])



figure
for i=1:59
    subplot(5,12,i)
    plot(mean(filters(:,i,:)/5e6,3))
    title(int2str(i))
    axis([1 500 -0.6 0.2])
end
subplot(8,8,61)
plot(mean(filters(:,30,:)/5e6,3))
axis([1 500 -1 0.5])
