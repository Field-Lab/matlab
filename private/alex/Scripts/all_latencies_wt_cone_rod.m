clear
load('/Users/alexth/Dropbox/paper/2_lightAdaptation/data/collect_param.mat','cone')

% select only units with parameter 9 (experiment type) equal to sel
% controlThr=5; % NO CONTROL CHECK
par=9;
sel=12;
col='rbk';
figure
set(gcf,'position',[52 522 1154 428])
subplot(2,1,1)
ms=10;
hold on
set(gca,'XTick',10:20:160,'xticklabel',{'8','7','6','5','4','3','2','1'})
xlabel('ND')
ylabel('latency, ms')

pulledLAT=cone(cone(:,49,9)==sel,:,2);

pulledAMP=abs(cone(cone(:,49,9)==sel,:,1));
all_ampl_nd4=pulledAMP(:,49:60);
sum(isnan(all_ampl_nd4)); % check for nonexisting values
meanND4=nanmean(all_ampl_nd4,2); % nanmean ampl at ND4 for each unit
sum(meanND4==0);

for i=1:size(pulledAMP,2)
    pulledAMP(:,i)=abs(pulledAMP(:,i)./meanND4);
end


j=1;
for i=1:20:160
    plot(i:i+11,pulledLAT(:,j:j+11),'.')  
    line([i+16,i+16],[0,200],'color','k')
    thisND=pulledLAT(:,j:j+11);
    thisND_AMP=pulledAMP(:,j:j+11);
    for k=1:12
        meanLat(j+k-1)=nanmean(thisND(thisND(:,k)>40&thisND(:,k)<160,k));        
        entityLat(j+k-1)=numel(thisND(thisND(:,k)>40&thisND(:,k)<160,k));
        stdLat(j+k-1)=nanstd(thisND(thisND(:,k)>40&thisND(:,k)<160,k))./sqrt(entityLat(j+k-1));

        meanAmp(j+k-1)=nanmean(thisND_AMP(:,k));
        entityAmp(j+k-1)=numel(thisND_AMP(:,k));
        stdAmp(j+k-1)=nanstd(thisND_AMP(:,k))./sqrt(entityAmp(j+k-1));

    end    
    j=j+12;
end

axis([0,160,0, 200])
title('cone only, Latency, experiments with 12 trials, n<=121')

j=1;
subplot(2,1,2)
hold on
for i=1:20:160
    errorbar(i:i+11,meanLat(j:j+11),stdLat(j:j+11),'color','r','linewidth',2)
    bar(i:i+11,entityLat(j:j+11),'facecolor',[0 0.7 1])
    line([i+16,i+16],[0,200],'color','k')  
    j=j+12;
end
axis([0,160,0, 200])

figure(5)
subplot(2,1,1)
j=1;
hold on
for i=1:20:160
    errorbar(i:i+11,meanLat(j:j+11),stdLat(j:j+11),'color','r','linewidth',2)
    line([i+16,i+16],[0,200],'color','k')  
    j=j+12;
end
axis([0,160,0, 200])
title('Latency; red - cones, blue - rods, black - wt')
set(gca,'XTick',10:20:160,'xticklabel',{'8','7','6','5','4','3','2','1'})
xlabel('ND')

figure(5)
subplot(2,1,2)
j=1;
hold on
for i=1:20:160
    errorbar(i:i+11,meanAmp(j:j+11),stdAmp(j:j+11),'color','r','linewidth',2)
    line([i+16,i+16],[0,200],'color','k')  
    j=j+12;
end
title('amplitude; red - cones, blue - rods, black - wt')
set(gca,'XTick',10:20:160,'xticklabel',{'8','7','6','5','4','3','2','1'})
xlabel('ND')



figure(6)
j=1;cnt=1;
for i=1:12:96
    subplot(2,4,cnt)
    hold on
    for k=0:11
        plot(pulledLAT(:,j+k),pulledAMP(:,j+k),'+','color',[1-k/14 0 0],'markersize',2)
    end
    axis([0,200,0, 3.5])
    j=j+12;
    cnt=cnt+1;
end




%%%%%%%% ROD ONLY %%%%%%%%%


clear
load('/Users/alexth/Dropbox/paper/2_lightAdaptation/data/collect_param.mat','rod')

% select only units with parameter 9 (experiment type) equal to sel
% controlThr=5; % NO CONTROL CHECK
par=9;
sel=12;
col='rbk';
figure
set(gcf,'position',[52 522 1154 428])
subplot(2,1,1)
ms=10;
hold on
set(gca,'XTick',10:20:160,'xticklabel',{'8','7','6','5','4','3','2','1'})
xlabel('ND')
ylabel('latency, ms')

pulledLAT=rod(rod(:,49,9)==sel,:,2);


pulledAMP=abs(rod(rod(:,49,9)==sel,:,1));
all_ampl_nd6=pulledAMP(:,25:36);
sum(isnan(all_ampl_nd6)); % check for nonexisting values
meanND6=nanmean(all_ampl_nd6,2); % nanmean ampl at ND6 for each unit
sum(meanND6==0);

for i=1:size(pulledAMP,2)
    pulledAMP(:,i)=abs(pulledAMP(:,i)./meanND6);
end


j=1;
for i=1:20:160
    plot(i:i+11,pulledLAT(:,j:j+11),'.')  
    line([i+16,i+16],[0,200],'color','k')
    thisND=pulledLAT(:,j:j+11);
    thisND_AMP=pulledAMP(:,j:j+11);
    for k=1:12
        meanLat(j+k-1)=nanmean(thisND(thisND(:,k)>40&thisND(:,k)<160,k));        
        entityLat(j+k-1)=numel(thisND(thisND(:,k)>40&thisND(:,k)<160,k));
        stdLat(j+k-1)=nanstd(thisND(thisND(:,k)>40&thisND(:,k)<160,k))./sqrt(entityLat(j+k-1));

        meanAmp(j+k-1)=nanmean(thisND_AMP(:,k));
        entityAmp(j+k-1)=numel(thisND_AMP(:,k));
        stdAmp(j+k-1)=nanstd(thisND_AMP(:,k))./sqrt(entityAmp(j+k-1));

    end    
    j=j+12;
end

axis([0,160,0, 200])
title('rod only, Latency, experiments with 12 trials, n<=53')

j=1;
subplot(2,1,2)
hold on
for i=1:20:160
    errorbar(i:i+11,meanLat(j:j+11),stdLat(j:j+11),'color','r','linewidth',2)
    bar(i:i+11,entityLat(j:j+11),'facecolor',[0 0.7 1])
    line([i+16,i+16],[0,200],'color','k')  
    j=j+12;
end
axis([0,160,0, 200])


figure(5)
subplot(2,1,1)
j=1;
hold on
for i=1:20:160
    errorbar(i:i+11,meanLat(j:j+11),stdLat(j:j+11),'color','b','linewidth',2)
    j=j+12;
end
axis([0,160,0, 200])




figure(5)
subplot(2,1,2)
j=1;
hold on
for i=1:20:160
    errorbar(i:i+11,meanAmp(j:j+11),stdAmp(j:j+11),'color','b','linewidth',2)
    line([i+16,i+16],[0,200],'color','k')  
    j=j+12;
end



figure(6)
j=1;cnt=1;
for i=1:12:96
    subplot(2,4,cnt)
    hold on
    for k=0:11
        plot(pulledLAT(:,j+k),pulledAMP(:,j+k),'o','color',[0 0 1-k/14],'markersize',2)
    end
    axis([0,200,0, 3.5])
    j=j+12;
    cnt=cnt+1;
end




%%%%%%%% WILD TYPE %%%%%%%%%


clear
load('/Users/alexth/Dropbox/paper/2_lightAdaptation/data/collect_param.mat','wt')

% select only units with parameter 9 (experiment type) equal to sel
% controlThr=5; % NO CONTROL CHECK
sel=12;
col='rbk';
figure
set(gcf,'position',[52 522 1154 428])
subplot(2,1,1)
ms=10;
hold on
set(gca,'XTick',10:20:160,'xticklabel',{'8','7','6','5','4','3','2','1'})
xlabel('ND')
ylabel('latency, ms')

pulledLAT=wt(wt(:,49,9)==sel,:,2);

pulledAMP=abs(wt(wt(:,49,9)==sel,:,1));
all_ampl_nd4=pulledAMP(:,49:60);
sum(isnan(all_ampl_nd4)); % check for nonexisting values
meanND4=nanmean(all_ampl_nd4,2); % nanmean ampl at ND4 for each unit
sum(meanND4==0);
for i=1:size(pulledAMP,2)
    pulledAMP(:,i)=abs(pulledAMP(:,i)./meanND4);
end



j=1;
for i=1:20:160
    plot(i:i+11,pulledLAT(:,j:j+11),'.')  
    line([i+16,i+16],[0,200],'color','k')
    thisND=pulledLAT(:,j:j+11);
    thisND_AMP=pulledAMP(:,j:j+11);
    for k=1:12
        meanLat(j+k-1)=nanmean(thisND(thisND(:,k)>40&thisND(:,k)<160,k));        
        entityLat(j+k-1)=numel(thisND(thisND(:,k)>40&thisND(:,k)<160,k));
        stdLat(j+k-1)=nanstd(thisND(thisND(:,k)>40&thisND(:,k)<160,k))./sqrt(entityLat(j+k-1));

        meanAmp(j+k-1)=nanmean(thisND_AMP(:,k));
        entityAmp(j+k-1)=numel(thisND_AMP(:,k));
        stdAmp(j+k-1)=nanstd(thisND_AMP(:,k))./sqrt(entityAmp(j+k-1));

    end     
    j=j+12;
end

axis([0,160,0, 200])
title('wt, Latency, experiments with 12 trials, n<=35')

j=1;
subplot(2,1,2)
hold on
for i=1:20:160
    errorbar(i:i+11,meanLat(j:j+11),stdLat(j:j+11),'color','r','linewidth',2)
    bar(i:i+11,entityLat(j:j+11),'facecolor',[0 0.7 1])
    line([i+16,i+16],[0,200],'color','k')  
    j=j+12;
end
axis([0,160,0, 200])



figure(5)
subplot(2,1,1)
j=1;
hold on
for i=1:20:160
    errorbar(i:i+11,meanLat(j:j+11),stdLat(j:j+11),'color','k','linewidth',2)
    j=j+12;
end
axis([0,160,0, 200])




figure(5)
subplot(2,1,2)
j=1;
hold on
for i=1:20:160
    errorbar(i:i+11,meanAmp(j:j+11),stdAmp(j:j+11),'color','k','linewidth',2)
    line([i+16,i+16],[0,200],'color','k')  
    j=j+12;
end
axis([0,160,0, 2])


figure(6)
j=1;cnt=1;
for i=1:12:96
    subplot(2,4,cnt)
    hold on
    for k=0:11
        plot(pulledLAT(:,j+k),pulledAMP(:,j+k),'x','color',[1-k/14 1-k/14 1-k/14],'markersize',2)
    end
    axis([0,200,0, 3.5])
    j=j+12;
    cnt=cnt+1;
end