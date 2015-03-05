2)	melanopsin lines and bars (mean +-sem), latencies of all trials at each ND
3)	melanopsin bars difference between first and last trial in each ND – latencies


pulled=mela(1:38,trialSelect,2);
clear
load('S:\Paper factory\Alex\2_lightAdaptation\data\collect_param.mat','mela')
trialSelect=sort([14:12:96 16:12:96 22:12:96 24:12:96]);
lat=mela(1:38,trialSelect,2);
ampl=mela(1:38,trialSelect,1);
polar=sign(sum(mela(1:38,trialSelect,4),2));
totTime=210;
timing=sort([2:30:totTime 7:30:totTime 23:30:totTime 29:30:totTime]);

subplot(2,1,1)
for i=1:28
    tmp=lat(:,i)<160&lat(:,i)>0;
    allMean(i)=mean(lat(tmp,i));
    allSem(i)=std(lat(tmp,i))/sqrt(sum(tmp));
end
errorbar(timing,allMean,allSem,'color','r')
set(gca,'XTick',15:30:totTime,'xticklabel',{'7','6','5','4','3','2','1'})
xlabel('ND')
for i=30:30:totTime-30
    line([i,i],[0,200],'color','k')
end
axis([0 totTime+1 80 160])
title('Melanopsin Latency, 3 retinas, 38 units, ND dur 30min')


subplot(2,1,2)
ampl=mela(1:38,trialSelect,1);
ampl=abs(ampl);
for i=1:38
    ampl(i,:)=ampl(i,:)/mean(ampl(i,5:8));  
%     ampl(i,:)=ampl(i,:)/max(ampl(i,:));
end
ampl(isinf(ampl))=NaN;
for i=1:28
    tmp=lat(:,i)<160&lat(:,i)>0;
    allMean(i)=nanmean(ampl(tmp,i));
    allSem(i)=nanstd(ampl(tmp,i))/sqrt(sum(tmp&~isnan(ampl(:,i))));
end
errorbar(timing,allMean,allSem,'color','r')
set(gca,'XTick',15:30:totTime,'xticklabel',{'7','6','5','4','3','2','1'})
xlabel('ND')
for i=30:30:totTime-30
    line([i,i],[0,1.3],'color','k')
end
axis([0 totTime+1 0 1.3])
title('Melanopsin Amplitude normalized to mean of nd6, 3 retinas, 38 units, ND dur 30min')







clear k n m
for i=1:size(pulled,2)
    pos=~isnan(polar(:,i));%&controlCheck(:,i);
    a=abs(pulled(pos,i));    
    a(isnan(a))=[];
    a(isinf(a))=[];
    k(i,1)=mean(a);
    k(i,2)=std(a)/sqrt(numel(a));
    m(i)=numel(a);
end
errorbar(k(:,1),k(:,2))

figure
for i=1:38
    plot(pulled(i,:),'.')    
    hold on
end
line([0,25],[0,0],'color','k')

for i=1:28
    l(i)=sum(pulled(:,i)>0&pulled(:,i)<200);
    k(i)=mean(pulled(pulled(:,i)>0&pulled(:,i)<200,i));
    m(i)=std(pulled(pulled(:,i)>0&pulled(:,i)<200,i))/sqrt(l(i));
end
a=sort([2:27:189 6:27:189 20:27:189 24:27:189]);
figure
errorbar(a,k,m)

figure
subplot(2,1,1)
hold on
for i=1:38
    plot(pulled(i,4:4:28))    
    plot(pulled(i,4:4:28),'.')        
end
line([0,8],[0,0],'color','k')
set(gca,'XTick',1:7,'xticklabel',{'7','6','5','4','3','2','1'})
subplot(2,1,2)
hold on
for i=1:38
    plot(pulled(i,1:4:28))    
    plot(pulled(i,1:4:28),'.')        
end
line([0,8],[0,0],'color','k')
set(gca,'XTick',1:7,'xticklabel',{'7','6','5','4','3','2','1'})



plot(k)

for i=4.5:4:28
    line([i,i],[0,200],'color','k')
end
clear a
l=1:4:28;
m=4:4:28;
a=abs(pulled(:,m))-abs(pulled(:,l));
p=a>-20&a<60;
for i=1:7
    plot(i,a(:,i),'.')
    line([0,7],[0,0],'color','k')
    hold on
    b(i)=mean(a(a(:,i)>-20&a(:,i)<60,i));
    k(i)=std(a(a(:,i)>-20&a(:,i)<60,i));
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


axis([12.5,96.5,0,1.5])
title('mela only, Amplitude, normalized to ND6, first trial, 2 experiments with 12 trials, n=~53, comb.On and OFF')
