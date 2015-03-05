
clear
dateR='20130629'; % date of the experiment, yyyymmdd

path2readSpikes=['/gpfs01/muench/data/Hartwig/MEA_data/',dateR,'/Spikes/'];
path2saveWF=['/gpfs01/muench/data/Hartwig/MEA_data/',dateR,'/waveforms/'];

if ~exist(path2saveWF,'file')
    mkdir(path2saveWF);
end

a=zeros(24,4000000,60);
s=dir([path2readSpikes,'*.mat']);
m=cell(length(s),60);
r=cell(length(s),60);

l=ones(60,1);
for i=1:length(s)
    i
    load([path2readSpikes,s(i).name])
    for j=1:60
        k=spikes{j};
        a(:,l(j):l(j)+size(k,2)-1,j)=k;
        m{i,j}=spikeTimes{j};
        r{i,j}=thr(j);
        l(j)=l(j)+size(k,2);
    end    
end

for i=1:60
    if i~=15
        spikes=a(:,1:l(i)-1,i);
        spikeTimes=m(:,i);
        thrs=r(:,i)';
        save([path2saveWF,'CH',int2str(i),'_',dateR],'spikes')
        save([path2saveWF,'CH',int2str(i),'_',dateR,'_spikeTimes'],'spikeTimes','thrs')
    end
    
end