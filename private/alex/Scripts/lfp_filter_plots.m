cd('/mnt/muench_data/user/alexandra/scripts')
date='20121017';

path2data=['/mnt/muench_data/data/alexandra/MEA_data/',date,'/LFP_filters/'];
keyWord='HCseq';
lfpFilters=dir([path2data,'*',keyWord,'*.mat']);

figure
subplot(2,1,2)
plot(all_ch(:,37:48)/59,'LineWidth',2)



figure
all_ch=0;
for ch=1:59
    load([path2data,lfpFilters(ch).name])
    all_ch=all_ch+HighFilter;
    subplot(6,10,ch)
    plot(HighFilter(:,[37,40,48]+12),'LineWidth',2)
end
subplot(6,10,60)
plot(all_ch(:,[37,40,48]+12)/59,'LineWidth',2)

nds='7654321234567654321'
figure
cnt=1;
for i=1:12:228
    subplot(4,5,cnt)
    plot(all_ch(1:500,i:i+11)/59,'LineWidth',2)
    axis([0 500 -30 20])
title(['ND',nds(cnt)])
cnt=cnt+1;
end




nds='654343212345'
figure
cnt=1;
for i=1:12:121
    subplot(3,4,cnt)
    if cnt<4
        plot(all_ch(1:500,i:i+11)/59,'LineWidth',2)
    elseif cnt==4
        plot(all_ch(1:500,i)/59,'LineWidth',2)
    else
        plot(all_ch(1:500,i-11:i)/59,'LineWidth',2)
    end
%     axis([0 500 -1000 1500])
title(['ND',nds(cnt)])
cnt=cnt+1;
end



t=[12 30:9:71 71 72:12:194];
tot_high=[1 25 9 9 9 9 9 1 12 12 12 12 12 9 9 9 9 9 9 9];
a=cumsum(tot_high);
a(end)=a(end)-1;
nds='9876543432123456783'
figure
cnt=1;
for i=1:19
    subplot(4,5,cnt)
    plot(all_ch(1:500,a(i):a(i+1)-1)/59,'LineWidth',2)

%     axis([0 500 -50 20])
title(['ND',nds(cnt)])
cnt=cnt+1;
end





nds='876543212345'
figure
cnt=1;
for i=1:4:40
    subplot(3,4,cnt)
    plot(all_ch(1:500,i:i+3)/59,'LineWidth',2)
%     axis([0 500 -1000 1500])
title(['ND',nds(cnt)])
cnt=cnt+1;
end


