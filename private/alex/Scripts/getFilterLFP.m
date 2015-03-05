clear
cd('/mnt/muench_data/user/alexandra/scripts')
date='20121017';

length2take=60; % time in s

path2data=['/mnt/muench_data/data/alexandra/MEA_data/',date,'/LFP/'];
keyWord='HCseq';
lfpData=dir([path2data,'*',keyWord,'*.mat']);

hekapath=['/mnt/muench_data/data/alexandra/MEA_data/',date,'/HEKA/'];
heka=dir([hekapath,'*.phys']);

file_list=[];
for i=1:length(heka)
    if ~isempty(regexp(heka(i).name,keyWord, 'once'))
        file_list=[file_list i];
    end
end
path2save=['/mnt/muench_data/data/alexandra/MEA_data/',date,'/LFP_filters/'];
if exist(path2save,'dir')    
    rmdir(path2save,'s')
end
mkdir(path2save);


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

for ch=1:59
    ch
    load([path2data,lfpData(ch).name])
    a=regexp(lfpData(ch).name,'CH');
    chs=lfpData(ch).name(a+2:a+3);
    if chs(2)=='_'
        chs=chs(1);
    end
    chs=str2num(chs);
    chs
    HighFilter=zeros(500,length(file_list));
    for trial_counter=1:length(file_list)
        %     trial_counter
        
        flips=protocol_accum(:,1,trial_counter); % flips in ms, rounded
        flips(flips>length2take*1000)=[];
        pxls=((protocol_accum(:,2,trial_counter)-30)/30); % pxs in ms, normalized -1:1
        pxls(length(flips)+1:end)=[];
        a=lfp(:,trial_counter);
        
        tmp=zeros(flips(end)+1,1); % 1 column - color, 2 column - spikes
        tmp([1; flips(1:end-1)])=pxls;
        
        % fill stimulus down to ms with brightness values where spikes
        % happen
        for i=1:50
            subscr=flips(1:end-1)+i;
            while length(tmp)<subscr(end)
                subscr(end)=[];
            end                
            subscr(tmp(subscr,1)~=0)=[];
            tmp(subscr,1)=tmp(subscr-1,1);
        end
        b=0;
        for i=1:length2take*1000-1000
            b=b+a(i+500)*tmp(i+499:-1:i);
        end
        b=b/(length2take*1000-1000);
        HighFilter(:,trial_counter)=b;   
    end
    save([path2save,'lfp_filters_CH',int2str(chs),'_',keyWord],'HighFilter')
end
