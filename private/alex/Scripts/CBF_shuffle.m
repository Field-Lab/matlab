cd('/mnt/muench_data/user/alexandra/scripts')
date='20120902_2';
unitspath=['/mnt/muench_data/user/alexandra/MEA_data/',date,'/easy_formatted_units/'];
units=dir([unitspath,'*CBFlicker_Shuffle.mat*']);

filter_length=500; % in ms


frames=ceil(filter_length/(100/6)); % in frames
load([unitspath,units(1).name]);
trials=size(spike_info.spike_times,1);

load(['/mnt/muench_data/user/alexandra/MEA_data/',date,'/checkers_sequences/seq1.mat'])
high_contr=(double(high_contr)-30)/30;

load(['/mnt/muench_data/user/alexandra/MEA_data/',date,'/checkers_sequences/coord_shift.mat'])


for units_cnt=1:length(units)
    units(units_cnt).name

    load([unitspath,units(units_cnt).name])
    
    SpikeCount=zeros(1,trials);
    CBFfilter=zeros(1600,filter_length,trials);
    tic
    for trial_counter=1:trials
        trial_counter
        tic
        % spike_info.name_info(trial_counter)
        flips=spike_info.flip_times{trial_counter}(:,1); % flips in ms, rounded
        pxls=spike_info.flip_times{trial_counter}(:,2); % pxs in ms, normalized -1:1
        spikes=spike_info.spike_times{trial_counter}'; % spikes in ms, rounded
        
        spikes(spikes<filter_length|spikes>flips(end-2))=[];
        spikes(spikes<flips(3)+filter_length)=[];
        length(spikes)
        SpikeCount(trial_counter)=length(spikes);
        
        starts_tmp=zeros(length(spikes),3);
        i=1;
        for cnt=spikes'
            beg=find(flips>cnt-filter_length,1)-1;
            take_fin=flips(beg+1)-(cnt-filter_length);            
            starts_tmp(i,:)=[pxls(beg) beg take_fin ];           
            i=i+1;
        end  
        for i=1:filter_length
            CBFfilter(:,i,trial_counter)=CBFfilter(:,i,trial_counter)+sum(high_contr(starts_tmp(:,1),:))'; % fill actual data
            starts_tmp(:,3)=starts_tmp(:,3)-1; % decrement "to take from this frame" value by 1
            starts_tmp(starts_tmp(:,3)==0,1)=starts_tmp(starts_tmp(:,3)==0,1)+1; % take the next frame for those frames which became 0
            tmp=flips(starts_tmp(starts_tmp(:,3)==0,2));
            starts_tmp(starts_tmp(:,3)==0,2)=starts_tmp(starts_tmp(:,3)==0,2)+1; % take the next flip for those frames which became 0
            starts_tmp(starts_tmp(:,3)==0,3)=flips(starts_tmp(starts_tmp(:,3)==0,2))-tmp; % take new "to take from this frame" values
        end     
        toc
    end
    toc
    save(['/mnt/muench_data/user/alexandra/MEA_data/',date,'/CBF_1000/',units(units_cnt).name(1:end-27),'_CBFlicker_LF.mat'],'CBFfilter','SpikeCount')
    
end





% prepare stimulus for trial 1

% get flip times (any unit)
load([unitspath,units(units_cnt).name])
flips=spike_info.flip_times{1}(:,1); % flips in ms, rounded
pxls=spike_info.flip_times{1}(:,2); % checker pattern number
spikes=spike_info.spike_times{1}';
spikes=spikes-flips(1);
spikes(spikes<=500)=[];
flips=flips-flips(1);
flips=flips(1:end-3);
spikes(spikes>flips(end))=[];

load(['/mnt/muench_data/user/alexandra/MEA_data/',date,'/checkers_sequences/seq1.mat'])
high_contr=(double(high_contr)-30)/30;

load(['/mnt/muench_data/user/alexandra/MEA_data/',date,'/checkers_sequences/coord_shift.mat'])
    

filter=zeros(42,42,500);
last_cnt=1;
dur=120000;
for snippets=1
    snippets
    frame=zeros(42,42,dur);
    cnt=last_cnt;
    start=1;
    while flips(cnt+1)<dur*snippets        
        coord=coord_shift(:,cnt);
        pxl=high_contr(cnt,:);
        pxl=reshape(pxl,13,13);
        tmp=zeros(840,840);
        topLeft=[30+coord(1),30+coord(2)];
        for i=1:60:13*60
            for j=1:60:13*60
                tmp(i+topLeft(1):i+topLeft(1)+59,j+topLeft(2):j+topLeft(2)+59)=pxl((i+59)/60,(j+59)/60);
            end
        end
        tmp1=zeros(42,42);
        % downsample by 3
        ii=1;
        for i=1:20:840
            jj=1;
            for j=1:20:840
                tmp1(ii,jj)=mean2(tmp(i:i+19,j:j+19));
                jj=jj+1;
            end
            ii=ii+1;
        end
        
        frame_length=flips(cnt+1)-flips(cnt);
        frame(:,:,start:start+frame_length-1)=repmat(tmp1,[1,1,frame_length]);
        start=start+frame_length;
        cnt=cnt+1;
        last_cnt=cnt;
    end
    
    a=find(spikes<dur*snippets&spikes>dur*(snippets-1)+500);
    
    for i=1:length(a)
        filter=filter+frame(:,:,spikes(a(i))-499-dur*(snippets-1):spikes(a(i))-dur*(snippets-1));
    end
end


imagesc(std(filter,0,3))
colormap gray


dat=reshape(filter(19,33,:),1,500);
plot(dat)
dat=reshape(filter,84*84,500);
plot(dat')

a=zeros(1,spikes(end));
for i=1:length(spikes)
    a(spikes(i)-500:spikes(i))=a(spikes(i)-500:spikes(i))+1;
end

b=zeros(1,flips(end));
for i=1:length(flips)-1
    b(flips(i)+1:flips(i+1))=pxls(i+1);
end


a=zeros(1,spikes(end));
a(spikes)=1;
a=reshape(a(1:605900),100,605900/100);
a=sum(a,1);
plot(a)
for i=1:length(spikes)
    a(spikes(i)-500:spikes(i))=a(spikes(i)-500:spikes(i))+1;
end