cd('/mnt/muench_data/user/alexandra/scripts')
date='20121023';

path2data=['/mnt/muench_data/user/alexandra/MEA_data/',date,'/LFP_filters/'];
keyWord='HCseq';
lfpFilters=dir([path2data,'*',keyWord,'*.mat']);

% path2save=['/mnt/muench_data/user/alexandra/MEA_data/',date,'/LFP_filters/'];
% if ~exist(path2save,'dir')    
%     mkdir(path2save);
% end

all_ch=0;
ind_comm=[];
for ch=1:59
    load([path2data,lfpFilters(ch).name])  
    [val,ind]=max(HighFilter(1:500,13:12:96));
    [val,ind1]=max(HighFilter(1:500,24:12:96));    
    kk(ch)=std(HighFilter(:,13));
    if kk(ch)>50
    ind_comm=[ind_comm; ind1-ind];
    end
    all_ch=all_ch+HighFilter;
end
save(['C:\Documents and Settings\atikidzhi\My Documents\Dropbox\paper\2_lightAdaptation\data\lfp_lat_diff_',date],'ind_comm')
figure
hold on
errorbar(mean(ind_comm),std(ind_comm)/sqrt(size(ind_comm,1)),'r')

nds='87654321'
figure
cnt=1;
for i=1:12:96
    subplot(3,4,cnt)
    plot(all_ch(1:500,i:i+11)/59,'LineWidth',2)
    axis([0 350 -1000 2500])
title(['ND',nds(cnt)])
cnt=cnt+1;
end



figure
all_ch=0;
for ch=1:59
    load([path2data,lfpFilters(ch).name])
    all_ch=all_ch+HighFilter;
    subplot(6,10,ch)
    plot(HighFilter(:,[13,8]),'LineWidth',2)
    title(int2str(kk(ch)))
end
subplot(6,10,60)
plot(all_ch(:,[13,8])/59,'LineWidth',2)



a=dir('C:\Documents and Settings\atikidzhi\My Documents\Dropbox\paper\2_lightAdaptation\data\lfp_lat_diff_*.mat');
c=cell(1,5);
all_ind=[];
for i=[1 3:5]
    load(a(i).name)
    if size(ind_comm,2)==7
        all_ind=[all_ind; ind_comm];
    elseif size(ind_comm,2)==4
        all_ind=[all_ind; [nan(49,3) ind_comm]];
    else
        all_ind=[all_ind; [nan(44,3) ind_comm nan(44,1)]];
    end
    c{i}=ind_comm;
end
all_mean=nanmean(all_ind);
all_std=nanstd(all_ind);
for i=1:7
    t(i)=sum(~isnan(all_ind(:,i)));
    b(i)=sqrt(t(i));
    all_std(i)=all_std(i)/b(i);
end
errorbar(all_mean,all_std,'r')





figure
set(gcf,'position',[1 31 1680 946])
bar(all_mean)
hold on
errorbar(all_mean,all_std,'.r')
for i=1:7
    text(i-0.3,all_mean(i)+all_std(i)+1,[int2str(all_mean(i)),'+-',num2str(all_std(i)),'ms  n=',int2str(t(i))]);
end
set(gca,'XTick',1:7,'xticklabel',{'7','6','5','4','3','2','1'})
axis([0 9 -5 50])


