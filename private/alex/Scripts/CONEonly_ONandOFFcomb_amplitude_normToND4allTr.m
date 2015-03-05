clear
load('S:\Paper factory\Alex\2_lightAdaptation\data\collect_param.mat','cone')

% select only units with parameter 9 (experiment type) equal to sel
% controlThr=5; % NO CONTROL CHECK
par=9;
sel=12;
col='rbk';
figure
set(gcf,'position',[1 31 1680 946])
ms=10;
hold on
set(gca,'XTick',6:12:96,'xticklabel',{'8','7','6','5','4','3','2','1'})

pulled=cone(:,:,1);
polar=cone(:,:,4);
all_ampl_nd4=pulled(:,49:60);
sum(isnan(all_ampl_nd4)); % check for nonexisting values
meanND4=nanmean(all_ampl_nd4,2); % mean ampl at ND6 for each unit
sum(meanND4==0);

clear k n m
for i=1:size(pulled,2)
    pos=~isnan(polar(:,i));%&controlCheck(:,i);
    a=abs(pulled(pos,i))./abs(meanND4(pos));
    k(i,1)=nanmean(a);
    k(i,2)=nanstd(a)/sqrt(sum(~isnan(a)));
    m(i,1)=numel(a);
end

a=find(~isnan(k(:,1)));
% k(isnan(k(:,1)),:)=[];



for i=13:12:96
    tt=a(a>=i&a<i+12);
    plot(tt,k(tt,1),'color',col(1),'linewidth',2)
    errorbar(tt,k(tt,1),k(tt,2),'.','color',col(2))
end

for i=12.5:12:96
    line([i,i],[0, 1.5],'color','k')
end
line([0,96.5],[0,0],'color','k')


axis([24.5,96.5,0, 1.5])
title('cone only, Amplitude, normalized to ND4 all trials,  experiments with 12/4 trials, n=~121, comb.On and OFF')
