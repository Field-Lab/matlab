clear
load('S:\Paper factory\Alex\2_lightAdaptation\data\collect_param.mat','wt')

% select only units with parameter 9 (experiment type) equal to sel
par=9;
sel=[4 6 12];

deltas=[2,12;1,11;1,12];
allExp=[];
for j=1:3
    unitSelect=sum(wt(:,:,par)==sel(j),2)>0;
    pulled=wt(unitSelect,:,2);
    allDeltas=[];
    for i=1:8
        allDeltas=[allDeltas pulled(:,deltas(j,2)+(i-1)*12)-pulled(:,deltas(j,1)+(i-1)*12)];
    end    
    allExp=[allExp; allDeltas];
end
clear k p
k(1:8,1)=nanmean(allExp);
k(1:8,2)=nanstd(allExp)./sqrt(sum(~isnan(allExp)));
for i=1:8
    [~,p(i)]=ttest(allExp(~isnan(allExp(:,i)),i));
end

figure
set(gcf,'position',[1 31 1680 946])
bar(k(:,1))
hold on
errorbar(k(:,1),k(:,2),'.r')
for i=1:8
    m(i)=sum(~isnan(allExp(:,i)));
    text(i-0.3,k(i,1)+k(i,2)+1,[int2str(k(i,1)),'+-',num2str(k(i,2)),'ms']);
    text(i-0.3,k(i,1)+k(i,2)+2.5,['p=',num2str(p(i))]);
    text(i-0.1,k(i,1)++k(i,2)+4,['n=',int2str(m(i))]);
end
set(gca,'XTick',1:8,'xticklabel',{'8','7','6','5','4','3','2','1'})
title('WT Latency differences between 1st and last trial, all experiments, p values and n for each bar')