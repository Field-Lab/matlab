function convert_wfs(waveform,spiketime)

global spikes spikeTimes

% make one array
long_spikes=zeros(34,10000000);
start=0;
for i=1:size(waveform,1)
    long_spikes(:,start+1:start+size(waveform{i},2))=waveform{i};
    start=start+size(waveform{i},2);
end
long_spikes(:,start+1:end)=[];

% alignment. minima at point 10.

spikes=zeros(24,start);
[~,pos]=min(long_spikes(11:19,:));
pos=pos+10;
for i=11:19
    first_point=i-9;
    last_point=i+14;
    spikes(:,pos==i)=long_spikes(first_point:last_point,pos==i);
end

spikeTimes=spiketime';
clear spiketime long_spikes