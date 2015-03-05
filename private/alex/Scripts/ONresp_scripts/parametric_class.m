clear
load('/Users/alexth/Desktop/old_stuff/ONresp/datarun_nd.mat')


for i=1:111
    tmp=datarun_nd(i).black.nd6;
    offBlack(i,1:3)=max(tmp(501:1000,:))-max(tmp(1:500,:));
end
a=offBlack(:,1)>5;
sum(a)
b=find(a);
totalCells=length(b);

for i=1:sum(a)
    subplot(8,8,i)
    plot(datarun_nd(b(i)).black.nd6)
end


figure
for i=1:sum(a)
    subplot(8,8,i)
    plot(datarun_nd(b(i)).lf.nd6)
end


%% ND7
% before 250
clear latMax latMin amMax amMin qc
for i=1:111
    tmp=datarun_nd(i).lf.nd7;
    [amMax(i),latMax(i)]=max(tmp(1:250,:));
    [amMin(i),latMin(i)]=min(tmp(1:250,:));
    qc(i)=std(tmp(1:250,:));    
end

minDif=abs(latMax-latMin);
thresh=5;
goods=minDif>30&qc>thresh;


figure
plot(latMax,latMin,'x')
hold on
plot(latMax(~goods(:)),latMin(~goods(:)),'+b')
xlabel('Lat Max, OFF')
ylabel('Lat Min, ON')
line([0 250],[0 250],'color','k','linewidth',2)
line([0 120],[120 120],'color','k')
line([120 120],[0 120],'color','k')


offs=(latMax-latMin>0)&goods&latMin>70;
ons=(-latMax+latMin>0)&goods&latMax>70;

figure(1)
subplot(2,4,1)
hold on
bar(mean(latMax(ons)));
errorbar(mean(latMax(ons)),std(latMax(ons)),'xr')
axis([0 4 0 200])
title(['ND7 ON cells, n=',int2str(sum(ons))])


subplot(2,4,5)
hold on
bar(mean(latMin(offs)));
errorbar(mean(latMin(offs)),std(latMin(offs)),'xr')
axis([0 4 0 200])
title(['ND7 OFF cells, n=',int2str(sum(offs))])



%% ND6
% before 250
clear latMax latMin amMax amMin qc
for i=1:111
    tmp=datarun_nd(i).lf.nd6;
    [amMax(i,1:3),latMax(i,1:3)]=max(tmp(1:250,:));
    [amMin(i,1:3),latMin(i,1:3)]=min(tmp(1:250,:));
    qc(i,1:3)=std(tmp(1:250,:));    
end

minDif=abs(latMax-latMin);
thresh=5;
goods=minDif>30&qc>thresh;


figure
plot(latMax,latMin,'x')
hold on
plot(latMax(~goods(:,1),1),latMin(~goods(:,1),1),'+b')
plot(latMax(~goods(:,2),2),latMin(~goods(:,2),2),'+g')
plot(latMax(~goods(:,3),3),latMin(~goods(:,3),3),'+r')
xlabel('Lat Max, OFF')
ylabel('Lat Min, ON')
line([0 250],[0 250],'color','k','linewidth',2)
line([0 120],[120 120],'color','k')
line([120 120],[0 120],'color','k')


offs=(latMax-latMin>0)&goods&latMin>70;
ons=(-latMax+latMin>0)&goods&latMax>70;

cond=ons(:,1)&ons(:,3);
[h,p]=ttest(latMax(cond,1),latMax(cond,3));
figure(1)
subplot(2,4,2)
hold on
bar(mean(latMax(cond,:)));
errorbar(mean(latMax(cond,:)),std(latMax(cond,:)),'xr')
if p<0.05
    k=max(mean(latMax(cond,:))+std(latMax(cond,:)))+10;
    line([1,1],[k,k+5],'color','k')
    line([3,3],[k,k+5],'color','k')
    line([1,3],[k+5,k+5],'color','k')
    if p<0.001
        text(2,k+5,'**')
    else
        text(2,k+5,'*')
    end
end
axis([0 4 0 200])
title(['ND6 ON cells, n=',int2str(sum(cond)),'   (c & wo, APB=',int2str(sum(ons(:,2))),')'])

cond=offs(:,1)&offs(:,2);
[h,p]=ttest(latMin(cond,1),latMin(cond,2));
cond=offs(:,1)&offs(:,3);
[h,p1]=ttest(latMin(cond,1),latMin(cond,3));
cond=offs(:,1)&offs(:,2)&offs(:,3);
subplot(2,4,6)
hold on
bar(mean(latMin(cond,:)));
errorbar(mean(latMin(cond,:)),std(latMin(cond,:)),'xr')
if p<0.05
    k=max(mean(latMin(cond,:))+std(latMin(cond,:)))+10;
    line([1,1],[k,k+5],'color','k')
    line([2,2],[k,k+5],'color','k')
    line([1,2],[k+5,k+5],'color','k')
    if p<0.001
        text(1.5,k+5,'**')
    else
        text(1.5,k+5,'*')
    end
end
if p1<0.05
    k=max(mean(latMin(cond,:))+std(latMin(cond,:)))+20;
    line([1,1],[k,k+5],'color','k')
    line([3,3],[k,k+5],'color','k')
    line([1,3],[k+5,k+5],'color','k')
    if p<0.001
        text(2,k+5,'**')
    else
        text(2,k+5,'*')
    end
end
axis([0 4 0 200])
title(['ND6 OFF cells, n=',int2str(sum(cond)),'   (c & APB & wo)'])



%% ND5
% before 250
clear latMax latMin amMax amMin qc
for i=1:111
    tmp=datarun_nd(i).lf.nd5;
    [amMax(i,1:3),latMax(i,1:3)]=max(tmp(1:250,:));
    [amMin(i,1:3),latMin(i,1:3)]=min(tmp(1:250,:));
    qc(i,1:3)=std(tmp(1:250,:));    
end

minDif=abs(latMax-latMin);
thresh=5;
goods=minDif>30&qc>thresh;


figure
plot(latMax,latMin,'x')
hold on
plot(latMax(~goods(:,1),1),latMin(~goods(:,1),1),'+b')
plot(latMax(~goods(:,2),2),latMin(~goods(:,2),2),'+g')
plot(latMax(~goods(:,3),3),latMin(~goods(:,3),3),'+r')
xlabel('Lat Max, OFF')
ylabel('Lat Min, ON')
line([0 250],[0 250],'color','k','linewidth',2)
line([0 100],[100 100],'color','k')
line([100 100],[0 100],'color','k')


offs=(latMax-latMin>0)&goods&latMin>70;
ons=(-latMax+latMin>0)&goods&latMax>70;

cond=ons(:,1)&ons(:,3);
[h,p]=ttest(latMax(cond,1),latMax(cond,3));
figure(1)
subplot(2,4,3)
hold on
bar(mean(latMax(cond,:)));
errorbar(mean(latMax(cond,:)),std(latMax(cond,:)),'xr')
if p<0.05
    k=max(mean(latMax(cond,:))+std(latMax(cond,:)))+10;
    line([1,1],[k,k+5],'color','k')
    line([3,3],[k,k+5],'color','k')
    line([1,3],[k+5,k+5],'color','k')
    if p<0.001
        text(2,k+5,'**')
    else
        text(2,k+5,'*')
    end
end
title(['ND5 ON cells, n=',int2str(sum(cond)),'   (c & wo, APB=',int2str(sum(ons(:,2))),')'])
axis([0 4 0 200])

cond=offs(:,1)&offs(:,2);
[h,p]=ttest(latMin(cond,1),latMin(cond,2));
cond=offs(:,1)&offs(:,3);
[h,p1]=ttest(latMin(cond,1),latMin(cond,3));
cond=offs(:,1)&offs(:,2)&offs(:,3);
subplot(2,4,7)
hold on
bar(mean(latMin(cond,:)));
errorbar(mean(latMin(cond,:)),std(latMin(cond,:)),'xr')
if p<0.05
    k=max(mean(latMin(cond,:))+std(latMin(cond,:)))+10;
    line([1,1],[k,k+5],'color','k')
    line([2,2],[k,k+5],'color','k')
    line([1,2],[k+5,k+5],'color','k')
    if p<0.001
        text(1.5,k+5,'**')
    else
        text(1.5,k+5,'*')
    end
end
if p1<0.05
    k=max(mean(latMin(cond,:))+std(latMin(cond,:)))+20;
    line([1,1],[k,k+5],'color','k')
    line([3,3],[k,k+5],'color','k')
    line([1,3],[k+5,k+5],'color','k')
    if p<0.001
        text(2,k+5,'**')
    else
        text(2,k+5,'*')
    end
end
axis([0 4 0 200])
title(['ND5 OFF cells, n=',int2str(sum(cond)),'   (c & APB & wo)'])


%% ND4
% before 250
clear latMax latMin amMax amMin qc
for i=1:111
    tmp=datarun_nd(i).lf.nd4;
    [amMax(i,1:3),latMax(i,1:3)]=max(tmp(1:250,:));
    [amMin(i,1:3),latMin(i,1:3)]=min(tmp(1:250,:));
    qc(i,1:3)=std(tmp(1:250,:));    
end

minDif=abs(latMax-latMin);
thresh=5;
goods=minDif>30&qc>thresh;


figure
plot(latMax,latMin,'x')
hold on
plot(latMax(~goods(:,1),1),latMin(~goods(:,1),1),'+b')
plot(latMax(~goods(:,2),2),latMin(~goods(:,2),2),'+g')
plot(latMax(~goods(:,3),3),latMin(~goods(:,3),3),'+r')
xlabel('Lat Max, OFF')
ylabel('Lat Min, ON')
line([0 250],[0 250],'color','k','linewidth',2)
line([0 100],[100 100],'color','k')
line([100 100],[0 100],'color','k')


offs=(latMax-latMin>0)&goods&latMin>70;
ons=(-latMax+latMin>0)&goods&latMax>70;

cond=ons(:,1)&ons(:,3);
[h,p]=ttest(latMax(cond,1),latMax(cond,3));
figure(1)
subplot(2,4,4)
hold on
bar(mean(latMax(cond,:)));
errorbar(mean(latMax(cond,:)),std(latMax(cond,:)),'xr')
if p<0.05
    k=max(mean(latMax(cond,:))+std(latMax(cond,:)))+10;
    line([1,1],[k,k+5],'color','k')
    line([3,3],[k,k+5],'color','k')
    line([1,3],[k+5,k+5],'color','k')
    if p<0.001
        text(2,k+5,'**')
    else
        text(2,k+5,'*')
    end
end
axis([0 4 0 200])
title(['ND4 ON cells, n=',int2str(sum(cond)),'   (c & wo, APB=',int2str(sum(ons(:,2))),')'])

cond=offs(:,1)&offs(:,2);
[h,p]=ttest(latMin(cond,1),latMin(cond,2));
cond=offs(:,1)&offs(:,3);
[h,p1]=ttest(latMin(cond,1),latMin(cond,3));
cond=offs(:,1)&offs(:,2)&offs(:,3);
subplot(2,4,8)
hold on
bar(mean(latMin(cond,:)));
errorbar(mean(latMin(cond,:)),std(latMin(cond,:)),'xr')
if p<0.05
    k=max(mean(latMin(cond,:))+std(latMin(cond,:)))+10;
    line([1,1],[k,k+5],'color','k')
    line([2,2],[k,k+5],'color','k')
    line([1,2],[k+5,k+5],'color','k')
    if p<0.001
        text(1.5,k+5,'**')
    else
        text(1.5,k+5,'*')
    end
end
if p1<0.05
    k=max(mean(latMin(cond,:))+std(latMin(cond,:)))+20;
    line([1,1],[k,k+5],'color','k')
    line([3,3],[k,k+5],'color','k')
    line([1,3],[k+5,k+5],'color','k')
    if p<0.001
        text(2,k+5,'**')
    else
        text(2,k+5,'*')
    end
end
axis([0 4 0 200])
title(['ND4 OFF cells, n=',int2str(sum(cond)),'   (c & APB & wo)'])
