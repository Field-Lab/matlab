clear
date='20120902_2';
path2filter=['S:\data\alexandra\MEA_data\',date,'\CBF_filter\'];
filters=dir([path2filter,'\*.mat']);
findHotSpot=1; % trial number

if ~exist(['S:\data\alexandra\MEA_data\',date,'\CBFlicker_LF\'],'file')
    mkdir(['S:\data\alexandra\MEA_data\',date,'\CBFlicker_LF\']);
end


for unit=1:length(filters)
    load([path2filter,filters(unit).name]);
    name=filters(unit).name(1:end-17);    
    
    tmp=CBFfilter(:,:,findHotSpot);
    tmp=tmp';
    [val co]=max(std(tmp));
    centrSize=3;
    surrSize=7;
    b=tmp(:,co);
    centr=[co-centrSize:co+centrSize];
    i=[];
    for j=-centrSize:centrSize
        i=[i centr+40*j];
    end
    i(i>1600)=[];
    i(i<1)=[];
    
    acc=zeros(size(CBFfilter,3),500);
    for trial=1:size(CBFfilter,3)
        tmp=CBFfilter(:,:,trial);
        tmp=tmp/max(abs(tmp(:)));
        acc(trial,1:500)=mean(tmp(i,:));
        acc(trial,1:500)=acc(trial,1:500)/max(abs(acc(trial,1:500)));
        acc(trial,1:500)=acc(trial,500:-1:1);        
    end
    HighFilters=acc';
    HighSpikeCount=SpikeCount;
    save(['S:\data\alexandra\MEA_data\',date,'\CBFlicker_LF\',name,'_filter'],'HighFilters','HighSpikeCount')
end