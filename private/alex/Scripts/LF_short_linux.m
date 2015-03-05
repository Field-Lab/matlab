clear
date='20130301';
codeWord='H30s'
unitspath=['/mnt/muench_data/data/alexandra/MEA_data/',date,'/easy_formatted_units/'];
units=dir([unitspath,'*',codeWord,'.mat']);
path2save=['/mnt/muench_data/data/alexandra/MEA_data/',date,'/LinearFilters/'];
if ~exist(path2save,'file')
    mkdir(path2save);
end
load([unitspath,units(1).name])
filter_length=500; %filter length in ms


for units_cnt=1:length(units)
    units(units_cnt).name
    if ~exist([path2save,units(units_cnt).name(1:end-27),'_',codeWord,'_linear_filter.mat'],'file')
        load([unitspath,units(units_cnt).name])
        trials=size(spike_info.flip_times,1);
        SpikeCount=zeros(1,length(trials));
        LinearFilter=zeros(filter_length,length(trials));

        for trial_counter=1:trials
            flips=spike_info.flip_times{trial_counter}(:,1); % flips in ms, rounded
            pxls=(((spike_info.flip_times{trial_counter}(:,2))-30)/30); % pxs in ms, normalized -1:1
            spikes=spike_info.spike_times{trial_counter}'; % spikes in ms, rounded
            spikes(spikes<filter_length|spikes>flips(end))=[];
            if ~isempty(spikes)
                tmp=zeros(flips(end)+1,2); % 1 column - color, 2 column - spikes
                tmp([1; flips(1:end-1)],1)=pxls;
                % fill stimulus down to ms with brightness values where spikes
                % happen
                for i=1:50
                    subscr=(flips(1:end-1)+i);
                    subscr(tmp(flips(1:end-1)+i,1)~=0)=[];
                    tmp(subscr,1)=tmp(subscr-1,1);
                end
                tmp(spikes,2)=1;
                n=zeros(length(spikes),filter_length); % filter
                for i=1:filter_length
                    n(:,i)=tmp(spikes-i+1,1);
                end
                n=[spikes, n];
 
                tmp=n(n(:,1)<flips(end-3)&n(:,1)>filter_length,2:end);
                    LinearFilter(1:filter_length,trial_counter)=sum(tmp);
                    SpikeCount(trial_counter)=size(tmp,1);
            end
        end
        save([path2save,units(units_cnt).name(1:end-27),'_FFFlicker_linear_filter.mat'],'HighFilters','HighSpikeCount','LowFilters','LowSpikeCount')
    end
end
