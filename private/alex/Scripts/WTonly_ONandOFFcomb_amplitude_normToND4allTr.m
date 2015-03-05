clear
load('/mnt/muench_data/Paper factory/Alex/2_lightAdaptation/data/collect_param.mat','wt')

pulled=wt(:,:,1);

a=nanmean(abs(pulled(:,12:12:96)));
cnt=1;
for i=12:12:96
    k(cnt)=sum(~isnan(pulled(:,i)));
    b(cnt)=nanstd(abs(pulled(:,i)))/sqrt(k(cnt));
    cnt=cnt+1;
end
errorbar(a,b)
set(gca,'XTick',1:8,'xticklabel',{'8','7','6','5','4','3','2','1'})
title('WT, last trial, amplitude mean+-sem')
xlabel('nd')
