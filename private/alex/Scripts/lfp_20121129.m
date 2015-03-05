cd('/mnt/muench_data/user/alexandra/scripts')
date='20121129';

path2data=['/mnt/muench_data/user/alexandra/MEA_data/',date,'/LFP_filters/'];
keyWord='HCseq';
lfpFilters=dir([path2data,'*',keyWord,'*.mat']);

% path2save=['/mnt/muench_data/user/alexandra/MEA_data/',date,'/LFP_filters/'];
% if ~exist(path2save,'dir')    
%     mkdir(path2save);
% end


figure
all_ch=0;
for ch=1:59
    load([path2data,lfpFilters(ch).name])
    all_ch=all_ch+HighFilter;
    subplot(6,10,ch)
    plot(HighFilter(:,[19,20]),'LineWidth',2)
    title(lfpFilters(ch).name(15:16))
end
subplot(6,10,60)
plot(all_ch(:,[19,20])/59,'LineWidth',2)

nds='765456754321357666'
figure
cnt=1;
for i=1:9:171
    subplot(4,5,cnt)
    plot(all_ch(1:500,i:i+8)/59,'LineWidth',2)
%     axis([0 500 -1000 1500])
title(['ND',nds(cnt)])
cnt=cnt+1;
end


nds='765456754321357666'
figure
cnt=1;
for i=1:9:171
    subplot(4,5,cnt)
    plot(HighFilter(1:500,i:i+8)/59,'LineWidth',2)
%     axis([0 500 -1000 1500])
title(['ND',nds(cnt)])
cnt=cnt+1;
end
