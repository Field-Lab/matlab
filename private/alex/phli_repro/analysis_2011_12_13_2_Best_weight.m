%% cell 223, off midget - compilation

mySingleConeCell=my_cells(my_cells(:,2)==5162,1); % cell ID 5301, cell index 165
spikes=datarun.spikes{165};
tr=datarun.triggers(1:2:end);

% contrast response

figure
colormap gray
for i=0:19
    if i<10
        map=load(['/Volumes/Analysis/2011-12-13-2/stimuli/2011-12-13-2_f08_allcones/map-000',int2str(i),'.txt']);
    else
        map=load(['/Volumes/Analysis/2011-12-13-2/stimuli/2011-12-13-2_f08_allcones/map-00',int2str(i),'.txt']);
    end
    subplot(4,5,i+1)
    imagesc(map+map_tmp)

end
myMaps=14:19; % starting from 0!

% calculate responses (# spikes to each contrast for each cone stimulation)
myMeanResponses=[];
for myCone=1:length(myMaps)
    myPresentations=find(stimulus.maps==myMaps(myCone));
    myContrasts=[];
    for i=1:length(myPresentations)
        myContrasts=[myContrasts stimulus.rgbs{myPresentations(i)}(:,1)];
    end
    
    for myCntr=1:4
        a=[];
        myCurrentContrast=find(myContrasts==allContrasts(myCntr));
        for i=1:length(myCurrentContrast)
            stimBegin=myPresentations(myCurrentContrast(i));
            tmp=tr(stimBegin);
            myResponse=spikes(spikes>=(tmp+0.05)&spikes<=(tmp+0.25))-tmp+0.05;
            a(i)=length(myResponse);
        end
        
        myMeanResponses(myCone,myCntr)=mean(a);
    end
end



% cone map with cone strength derived from relative response strength (NOT scaling)
myNormResponse=myMeanResponses/max(myMeanResponses);
maps_comb=0;
for i=1:length(myMaps)
    map=load(['/Volumes/Analysis/2011-12-13-2/stimuli/2011-12-13-2_f08_allcones/map-00',int2str(myMaps(i)),'.txt']);
    map(map>0)=myNormResponse(i);
    maps_comb=maps_comb+map;
end
figure
colormap gray
imagesc(-maps_comb(420:480,432:492))
title('Cell ID 223')


% RF input strength

%plot sta of this cell
sta=squeeze(datarunA.stas.stas{223});
mySTA=imresize(sta(:,:,5),2);
figure
colormap gray
imagesc(mySTA)

 % starting from 0!
myConeWeight=[];
myMap=zeros(61,61);
for myCone=1:length(myMaps)
    if myMaps(myCone)<10
        map=load(['/Volumes/Analysis/2011-12-13-2/stimuli/2011-12-13-2_f08_allcones/map-000',int2str(myMaps(myCone)),'.txt']);
    else
        map=load(['/Volumes/Analysis/2011-12-13-2/stimuli/2011-12-13-2_f08_allcones/map-00',int2str(myMaps(myCone)),'.txt']);
    end
    
    tmp=-mySTA.*map;
    tmp=tmp(420:480,432:492);
    myConeWeight(myCone)=sum(tmp(:));
    tmp(tmp~=0)=sum(tmp(:));
    myMap=myMap+tmp;

end

figure
colormap gray
imagesc(-myMap)

rfInputStrength=myConeWeight/max(myConeWeight);


% first calculate minimum squared error between RF cone inputs and simple
% relative # spikes for single cones to get initial estimates for fitting


% find strongest cone
strongestCone=find(myMeanResponses(:,1)==max(myMeanResponses(:)));
myNormResponse=myMeanResponses(:,1)/myMeanResponses(strongestCone,1);
coneInputStrength=myNormResponse;

myFun=fittype('a*x');
rfFit=fit(coneInputStrength,rfInputStrength',myFun);
totalScaleFactor=rfFit.a;

% calculate cone map
% initial scale factors
scales=coneInputStrength*totalScaleFactor;
myNormResponses=myMeanResponses/myMeanResponses(strongestCone,1);
myNormResponses=myNormResponses(:,end:-1:1);
myNormResponses=myNormResponses';
tmp=myNormResponses(:);

x=(0.25:0.25:1)';
xdata=[];
for k=1:size(myMeanResponses,1)
    xdata=[xdata; x*scales(k)];
end

gausCum=fittype('cdf(''norm'',x,mu,sigma)*height');
[firstFit,gof]=fit(xdata,tmp,gausCum,'Start',[1.3,0.7,0.5]);


for cnt=1:10
    mu=firstFit.mu;
    height=firstFit.height;
    sigma=firstFit.sigma;
    gausCumScales=eval(['fittype(''cdf(''''norm'''',scale*x,',num2str(mu),',',num2str(sigma),')*',num2str(height),''')']);
    for i=1:size(myMeanResponses,1)
        secondFit=fit((0.25:0.25:1)',myNormResponses(:,i),gausCumScales,'Start',scales(i));
        scales(i)=secondFit.scale;
    end
    xdata=[];
    for k=1:size(myMeanResponses,1)
        xdata=[xdata; x*scales(k)];
    end
    [firstFit,gof]=fit(xdata,tmp,gausCum,'Start',[height mu sigma]);
    gof.rsquare    
end




figure
set(gcf,'position',[1621 525 526 481])
hold on
x=0.25:0.25:1;
markers='+o^sxd...';
for i=1:size(myMeanResponses,1)
    plot(scales(i)*-x,myNormResponses(:,i)*max(myMeanResponses(:)),[markers(size(myMeanResponses,1)-i+1),'k'],'markersize',15,'linewidth',2)
end
axis tight
a=get(gca,'Xlim');
axis([a(1)*1.1 0 0 7])
hold on
x=-1:0.05:0;
plot(x,cdf('norm',-x,firstFit.mu,firstFit.sigma)*firstFit.height*max(myMeanResponses(:)),'k');

title({'contrast-response function, cell ID 5301 (223)',['SCALES: ', num2str(scales',3)]},'fontsize',18)
set(gca,'xtick',0,'xticklabel','0','fontsize',26)
set(gca,'ytick',0:2:6,'yticklabel',{'0','2','4','6'},'fontsize',26)
ylabel('# of spikes','fontsize',28)
xlabel('contrast (arbitrary units)','fontsize',28)





figure
set(gcf,'position',[1870         469         496         452])
plot(rfInputStrength,scales,'sk','markersize',15,'linewidth',2)
axis([0 1.1 0 1.1])
line([0,1.1],[0,1.1],'linestyle','--','color','k')

title('cell ID 223, Fig6G','fontsize',18)
set(gca,'xtick',0:0.2:1,'fontsize',26)
set(gca,'ytick',0:0.2:1,'fontsize',26)
ylabel('contrast response scaling','fontsize',24)
xlabel('receptive field input strength','fontsize',24)

%% Best weights


%plot sta of this cell
sta=squeeze(datarunA.stas.stas{223});
tmp=-sta(:,:,5);
[a,b]=find(tmp>robust_std(tmp(:))*3.5);
tc=0;
for i=1:length(a)
    tc=tc+sta(a(i),b(i),:);
end
tc=squeeze(tc)'/length(a);
figure
plot(tc)

tmp=reshape(sta,300*300,6);
tmp=tmp.*repmat(tc,300*300,1);
tmp=sum(tmp,2);
tmp=reshape(tmp,300,300);
figure
colormap gray
imagesc(-tmp)
resizedTC=imresize(-tmp,2, 'method','nearest');

figure
colormap gray
imagesc(resizedTC(420:480,432:492))
mySTA=resizedTC;


resizedSTA=zeros(600,600,6);
for i=1:6
    resizedSTA(:,:,i)=imresize(sta(:,:,i),2, 'method','nearest');
end

figure
colormap gray
imagesc(resizedSTA(:,:,5))

mySTA=resizedSTA(:,:,5);
figure
colormap gray
imagesc(mySTA(420:480,432:492))



 % starting from 0!
myConeWeight=[];
myMap=zeros(61,61);
figure
for myCone=1:length(myMaps)
    if myMaps(myCone)<10
        map=load(['/Volumes/Analysis/2011-12-13-2/stimuli/2011-12-13-2_f08_allcones/map-000',int2str(myMaps(myCone)),'.txt']);
    else
        map=load(['/Volumes/Analysis/2011-12-13-2/stimuli/2011-12-13-2_f08_allcones/map-00',int2str(myMaps(myCone)),'.txt']);
    end
    
    tmp=-mySTA.*map;
    tmp=tmp(420:480,432:492);
    
    
    subplot(2,3,myCone)
    colormap gray
    imagesc(-tmp)
    
    myConeWeight(myCone)=sum(tmp(:));
    tmp(tmp~=0)=sum(tmp(:));
    myMap=myMap+tmp;

end

figure
colormap gray
imagesc(-myMap)

scales=scales/max(scales);


rfInputStrength=myConeWeight/max(myConeWeight);


% all cones
figure
myTC=cell(6,1);
myTCsorted=myTC;myTC4sorted=myTC;
for myCone=1:6;
    map=load(['/Volumes/Analysis/2011-12-13-2/stimuli/2011-12-13-2_f08_allcones/map-00',int2str(myMaps(myCone)),'.txt']);

    tmp=-mySTA.*map;
    tmp=tmp(420:480,432:492);
    k=find(tmp);        
    for i=1:6
        tt=resizedSTA(420:480,432:492,i);
        myTC{myCone}=[myTC{myCone} tt(k)];        
    end
    
    subplot(2,3,myCone)
    plot(myTC{myCone}')
    title(['cone ',int2str(myCone)])
    
    for myTCpoint=1:6
        myTCsorted{myCone, myTCpoint}=-sort(myTC{myCone}(:,myTCpoint));
    end
    
    
end

a=round(rand(1500,1)*(300*300));
for j=1:6
    tt=resizedSTA(:,:,j);
    for i=1:length(a)
        myTCcontrol(i,j)=tt(a(i));
    end
end


figure
for myTCpoint=1:6
    myRatio=[];
    for i=1:74
        weight=[];
        for j=1:6
            weight(j)=abs(sum(myTCsorted{j,myTCpoint}(1:i)));%-sum(myTC4sorted{j,myTCpoint}(1:i));
            
        end
        %     weight=weight/(sum(weight));
        myweights=weight/max(weight);
        myRatio=[myRatio; myweights];
    end
    subplot(2,3,myTCpoint)
    plot(myRatio)
    hold on
    for i=1:6
        line([0,74],[scales(i), scales(i)],'color','k','linewidth',2,'linestyle','--')
        %     line([0,74],[rfInputStrength(i), rfInputStrength(i)],'color','m','linewidth',2)
    end
    axis([0 75 0 1.1])
    title(int2str(myTCpoint))
end



myFun=fittype('a*x');
rfFit=fit(coneInputStrength,rfInputStrength',myFun);
totalScaleFactor=rfFit.a;

figure
for myTCpoint=1:6
    bestPoint=[];
    for i=1:74
        for myCone=1:6
            weight(myCone)=abs(sum(myTCsorted{myCone,myTCpoint}(1:i)));%-sum(myTCsorted{j,myTCpoint}(1:i));
        end
        myweights=weight/max(weight);
        bestPoint(i)=sqrt(sum((scales-myweights').^2));
        
    end
    
    subplot(2,3,myTCpoint)
    plot(bestPoint,'-*')
    title(int2str(myTCpoint))
end



figure
for i=1:74
    for j=1:6
        weight(j)=sum(myTCsorted{j,4}(1:i))-sum(myTCsorted{j,3}(1:i));
    end
    myweights=weight/max(weight);
    bestPoint(i)=sqrt(sum((scales-myweights').^2));
    
end

plot(bestPoint(1:60),'-*')

figure
for i=1:6
    subplot(2,3,i)
    plot(myTCcontrol')
    hold on
    a=[];
    for j=1:6
        a=[a myTCsorted{i,j}];
    end
     plot(a','k','linewidth',2)
end


figure
for i=1:6
    subplot(2,3,i)
    a=myTCcontrol-repmat(mean(myTCcontrol(:,1:2),2),1,6);
    plot(a')
    hold on
    a=[];
    for j=1:6
        a=[a myTCsorted{i,j}];
    end
    b=a-repmat(mean(a(:,1:2),2),1,6);
    plot(b','k','linewidth',2)
end


figure
for i=1:6
    subplot(2,3,i)
    hold on
    a=[];
    for j=1:6
        a=[a myTCsorted{i,j}];
    end
    
    kk=a-repmat(mean(a(:,4:5),2),1,6);
%     b=b./abs(repmat(b(:,2),1,6));
    a=std(kk(:,2));  
    
    b=myTCcontrol-repmat(mean(myTCcontrol(:,4:5),2),1,6);
%     b=b./abs(repmat(b(:,2),1,6));
%     c=std(b(:,2));
%     b=b/(c/a);
    plot(b')
    
    plot(kk','k','linewidth',2)
    
%     axis([1 6 -1000 1000])
end


%%

[cell_list, failed_cells]=map_ei(datarunB,datarunA);

datarunA = load_data('/Volumes/Analysis/2011-12-13-2/streamed/data008-0/data008-0');
opt=struct('verbose',1,'load_params',1,'load_neurons',1,'load_obvius_sta_fits',false,'load_ei',true);
datarunA=load_data(datarunA,opt);
datarunA = load_params(datarunA,struct('verbose',1));


datarunB = load_data('/Volumes/Analysis/2011-12-13-2/streamed/data000-0/data000-0');
opt=struct('verbose',1,'load_params',1,'load_neurons',1,'load_obvius_sta_fits',false,'load_ei',true);
datarunB=load_data(datarunB,opt);
datarunB = load_params(datarunB,struct('verbose',1));



%% sigmoid


x=-5:0.1:5;
y=ones(size(x))./(1+exp(-x));
figure
plot(x,-y)

a=coneInputStrength/1.19

