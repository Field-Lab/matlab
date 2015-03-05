cd('/mnt/muench_data/user/alexandra/scripts')
date='20120902_2';
unitspath=['/mnt/muench_data/data/alexandra/MEA_data/',date,'/easy_formatted_units/'];
units=dir([unitspath,'*CBFlicker.*']);

filter_length=500; % in ms
path2save=['/mnt/muench_data/data/alexandra/MEA_data/',date,'/CBF_filter/'];

if ~exist(path2save,'dir')
    mkdir(path2save);
end

frames=ceil(filter_length/(100/6)); % in frames
load([unitspath,units(1).name]);
trials=size(spike_info.spike_times,1);

% load(['/mnt/muench_data/user/alexandra/MEA_data/',date,'/HC_seq.mat'])
load('/mnt/muench_data/user/alexandra/stimfiles/_hlc2_/HC_seq.mat')
high_contr=(double(high_contr)-30)/30;

for units_cnt=1:length(units)
    units(units_cnt).name
    if ~exist([path2save,units(units_cnt).name(1:end-27),'_CBFlicker_LF.mat'])
    load([unitspath,units(units_cnt).name])
    
    SpikeCount=zeros(1,trials);
    CBFfilter=zeros(1600,filter_length,trials);
    tic
    for trial_counter=1:trials
        trial_counter
        tic
        % trial_counter=28;
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
    save([path2save,units(units_cnt).name(1:end-27),'_CBFlicker_LF.mat'],'CBFfilter','SpikeCount')
    end
end