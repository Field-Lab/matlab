function get_spikes_v2(ch,date)
pathway=['/mnt/muench_data/data/alexandra/MEA_data/',date,'/MEA_bin/channels/channel',int2str(ch),'/'];
file_list=dir([pathway,'*.mat']);


if ~exist(['/mnt/muench_data/data/alexandra/MEA_data/',date,'/waveforms'],'file')
    mkdir(['/mnt/muench_data/data/alexandra/MEA_data/',date,'/waveforms']);
end


% HIGH PASS - for SPIKES
passH=500;%frequency of the high-pass filter
sf=25000;%sampling frequency
fNormH = passH/(sf/2);
orderH=10;%filter order
[aH,bH] = butter(orderH, fNormH, 'high');
clear sf passH fNormH orderH

maxSpikesSecond=150;


% allocate memory
spikeTimes=cell(1,length(file_list));
thrs=zeros(length(file_list),1);
spikes=zeros(24,6000000);
last=0; % last spike processed

for fileID=1:length(file_list)
    load([pathway, file_list(fileID).name]);
    
    %getting spikes - each cell for stimulus file; only time stamps
    %of spikes
    data=filtfilt(aH,bH,data);
    recDuration=length(data)/25000; % duration of recording in s
    maxSpikesAssumed=round(maxSpikesSecond*recDuration); % max spike amount to consider
    
    tmp=[0 diff(sign(diff(data)))];
    tmp=sort(data(tmp==2));
    thr=tmp(ceil(maxSpikesAssumed));
    tmp=diff(data>thr);
    timeStamps=find(tmp==-1)+1;
    
    % remove spikes very close to the beginning and the end
    timeStamps(timeStamps<20)=[];
    timeStamps(timeStamps>length(data)-30)=[];
    
    thrs(fileID)=thr;
    waveform=zeros(24,length(timeStamps));
    for k=-9:14
        waveform(k+10,:)=data(timeStamps+k);
    end
    
    timeStamps=timeStamps/25000; % spike times in s
    spikes(:,last+1:last+length(timeStamps))=waveform; % waveforms
    spikeTimes{fileID}=timeStamps;
    last=last+length(timeStamps);
 
end
spikes(:,last+1:end)=[];
last
save(['/mnt/muench_data/data/alexandra/MEA_data/',date,'/waveforms/CH',int2str(ch),'_',date,'.mat'],'spikes');
save(['/mnt/muench_data/data/alexandra/MEA_data/',date,'/waveforms/CH',int2str(ch),'_',date,'_spikeTimes.mat'],'spikeTimes','thrs')
