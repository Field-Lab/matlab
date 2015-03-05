clear
date='20130220_1';
unitspath=['/mnt/muench_data/data/alexandra/MEA_data/',date,'/easy_formatted_units/'];
units=dir([unitspath,'*FFFlicker.mat']);

path2save=['/mnt/muench_data/data/alexandra/MEA_data/',date,'/FFFlicker_LF_raw/'];

if ~exist(path2save,'file')
    mkdir(path2save);
end
load([unitspath,units(1).name])
filter_length=500; %filter length in ms
hseq=[];%m=[];
for i=1:length(spike_info.name_info)
    if ~isempty(cell2mat(regexp(spike_info.name_info(i),'HCseq')))
        hseq=[hseq i];
    end
end



cnt=1;
trial_counter=28;
flips=spike_info.flip_times{trial_counter}(:,1); % flips in ms, rounded
pxls=(((spike_info.flip_times{trial_counter}(:,2))-30)/30); % pxs in ms, normalized -1:1
spikes=spike_info.spike_times{trial_counter}'; % spikes in ms, rounded
spikes(spikes<flips(1)+filter_length|spikes>flips(end-2)-filter_length)=[];

HighSpikeCount(cnt)=length(spikes);

% coarseFilter=zeros(30,length(spikes));
% for i=1:length(spikes)
%     a=find(flips>spikes(i),1);
%     k(i)=a;
%     coarseFilter(1:30,i)=pxls(a:-1:a-29);
% end


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

rawStim=zeros(length(spikes),filter_length); % filter
for i=1:filter_length
    rawStim(:,i)=tmp(spikes-i+1,1);
end

% save([path2save,units(units_cnt).name(1:end-27),'_FFFlicker_linear_filter_raw_double_spikes.mat'],'HighFilterRow','HighSpikeCount','LowFilterRow','LowSpikeCount')


a=convolved(spikes,40,flips(end));
a=a(120:end-120);
tt=mean(rawStim);
plot(tt)


clear ind loc x

kk=tt(245:-1:1);
m=conv(tmp(:,1),kk,'same');
b=m.*tmp(:,1);

for i=-100:300
    x(i+101)=corr(b(6000:55000), a(6000+i:55000+i)');
end
figure
plot(-100:300,x)

[x,y]=xcorr(a(6000:55000)',b(6000:55000),'coeff');
figure
plot(y,x)
