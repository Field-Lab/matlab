clear
clear
dateR='20130302_2'; % date of the experiment, yyyymmdd

path2readSpikes=['/mnt/muench_data/data/alexandra/MEA_data/',dateR,'/Spikes/'];
path2saveWF=['/mnt/muench_data/data/alexandra/MEA_data/',dateR,'/waveforms/'];
path2savePCA=['/mnt/muench_data/data/alexandra/MEA_data/',dateR,'/pca/'];

if ~exist(path2saveWF,'file')
    mkdir(path2saveWF);
end
if ~exist(path2savePCA,'file')
    mkdir(path2savePCA);
end


a=zeros(24,4000000,60);
s=dir([path2readSpikes,'*.mat']);
m=cell(length(s),60);
r=cell(length(s),60);
pca1=cell(length(s),60);
pca2=cell(length(s),60);

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
        pca1{i,j}=pc1{j};
        pca2{i,j}=pc2{j};
    end    
end

for i=1:60
    if i~=15
        spikes=a(:,1:l(i)-1,i);
        spikeTimes=m(:,i);
        thrs=r(:,i)';
        pc1s=pca1(:,i);
        pc2s=pca2(:,i);
        save([path2saveWF,'CH',int2str(i),'_',dateR],'spikes')
        save([path2saveWF,'CH',int2str(i),'_',dateR,'_spikeTimes'],'spikeTimes','thrs')
        save([path2savePCA,'CH',int2str(i),'_',dateR],'pc1s','pc2s')
    end
    
end