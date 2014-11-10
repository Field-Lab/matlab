stimulus=read_stim_lisp_output_ath('2011-12-13-2','s10');

% parsed=parse_stim_rgbs_ath(stimulus);
% 
% myUniqueCone=[];
% for i=1:length(stimulus.rgbs)
%     
%     a=stimulus.rgbs{i}(:,1);
%     if nnz(a)==1 && a(1)~=0 
%         myUniqueCone=[myUniqueCone sum(a)];
%     end
%     
% end


map_all=0;
for i=1:19
    if i<10
        map=load(['/Volumes/Analysis/2011-12-13-2/stimuli/2011-12-13-2_f08_allcones/map-000',int2str(i),'.txt']);
    else
        map=load(['/Volumes/Analysis/2011-12-13-2/stimuli/2011-12-13-2_f08_allcones/map-00',int2str(i),'.txt']);
    end
    map_all=map_all+map*i;

end
figure
colormap gray
imagesc(map_all)

%% load datarun with single flashes

piece = '2011-12-13-2';
run = 'data010';

datarun = load_data(['/Volumes/Analysis/' piece '/' run '/' run]);
opt = struct('verbose',1,'load_params',1,'load_neurons',1,'load_obvius_sta_fits',false,'load_ei',true);
datarun = load_data(datarun,opt);
datarun = load_params(datarun,struct('verbose',1)); 


%% load datarun with white noise
piece = '2011-12-13-2';
run = 'data008-0';

% define data path
datarunA = load_data('/Volumes/Analysis/2011-12-13-2/streamed/data008-0/data008-0');
opt=struct('verbose',1,'load_params',1,'load_neurons',1,'load_obvius_sta_fits',false,'load_ei',true);
datarunA=load_data(datarunA,opt);
datarunA = load_params(datarunA,struct('verbose',1));
datarunA = load_sta(datarunA,'load_sta',[]);
datarunA = set_polarities(datarunA);
[datarunA, ~] = import_single_cone_data(datarunA, '/Volumes/Analysis/2011-12-13-2/_snle_acquisition_2011-12-13-2_data008-0_data008-0-bayes-msf_15.00-BW-2-6');
datarunA = make_mosaic_struct(datarunA);
datarunA = get_sta_fits_from_vision(datarunA);  
datarunA = make_voronoi_masks(datarunA);


%% map cells in dataruns (single flashes as master)

[cell_list, failed_cells]=map_ei(datarun,datarunA);

my_cells=zeros(length(cell_list),2);
for i=1:length(cell_list)
    my_cells(i,1)=datarun.cell_ids(i);
    if isempty(cell_list{i})
        my_cells(i,2)=0;
    else
        my_cells(i,2)=cell_list{i};
    end
end

%%

imagesc(map_all)
plot_rf_summaries(datarunA, {3}, 'clear', false, 'scale',2,'label', true, 'label_color', 'r', 'plot_fits', true, 'fit_color', 'r')

datarunA = load_sta(datarunA,'load_sta','all');

map_tmp=map_all;
map_tmp(map_tmp>0)=0.3;
for i=291:330
    sta=imresize(double(squeeze(datarunA.stas.stas{i}(:,:,:,5))),2);    
    if abs(min(sta(:)))>max(sta(:)) % off cell
        sta=-sta;
    end
    sta=sta/max(abs(sta(:)));
    sta(sta<0)=0;
    
    comb=zeros(600,600,3);
    comb(:,:,1)=sta;
    comb(:,:,2)=map_tmp;
    figure
    imagesc(comb)
    title([int2str(i),'  cell ID', int2str(datarunA.cell_ids(i))])
end

cellIndices=[296 297 220 205 17 223 161 132 94 95 38 55 67 61];

for i=cellIndices
    sta=imresize(double(squeeze(datarunA.stas.stas{i}(:,:,:,5))),2);    
    if abs(min(sta(:)))>max(sta(:)) % off cell
        sta=-sta;
    end
    sta=sta/max(abs(sta(:)));
    sta(sta<0)=0;
    
    comb=zeros(600,600,3);
    comb(:,:,1)=sta;
    comb(:,:,2)=map_tmp;
    figure
    imagesc(comb)
    title([int2str(i),'  cell ID', int2str(datarunA.cell_ids(i))])
end


subIndices=[223 161 61 67]; % 67 is additional cell, not shown on Fig 6F
cnt=1;
figure
for i=subIndices
    sta=imresize(double(squeeze(datarunA.stas.stas{i}(:,:,:,5))),2);    
    if abs(min(sta(:)))>max(sta(:)) % off cell
        sta=-sta;
    end
    sta=sta/max(abs(sta(:)));
    sta(sta<0)=0;
    
    comb=zeros(600,600,3);
    comb(:,:,1)=sta;
    comb(:,:,2)=map_tmp;
    subplot(2,2,cnt)
    imagesc(comb)
    title([int2str(i),'  cell ID', int2str(datarunA.cell_ids(i))])
    cnt=cnt+1;
end

allContrasts=unique(myContrasts);



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
mySTA=imresize(sta(:,:,5),2, 'method','nearest');
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

fig6G=imread('/Users/alexth/Desktop/Peters_data_replica/fig6G_cellB.jpeg');
figure
colormap gray
imagesc(fig6G)
coord=[507 232; 430 241; 327 325; 290 340; 234 444; 201 498];
x01=[142 507];
y01=[540 172];
rfs=(coord(:,1)-x01(1))/diff(x01);
rfs=rfs(end:-1:1);
rfs-rfInputStrength'
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
axis([a(1)*1.1 a(2)*1.1 0 7])

title({'contrast-response function, cell ID 5301, Norm',['SCALES: ', num2str(scales,3)]},'fontsize',18)
set(gca,'xtick',0,'xticklabel','0','fontsize',26)
set(gca,'ytick',0:2:6,'yticklabel',{'0','2','4','6'},'fontsize',26)
ylabel('# of spikes','fontsize',28)
xlabel('contrast (arbitrary units)','fontsize',28)





figure
set(gcf,'position',[1870         469         496         452])
plot(rfInputStrength,scales,'sk','markersize',15,'linewidth',2)
axis([0 1.1 0 1.1])
line([0,1.1],[0,1.1],'linestyle','--','color','k')

title('cell ID 161, Fig6G, norm fit','fontsize',18)
set(gca,'xtick',0:0.2:1,'fontsize',26)
set(gca,'ytick',0:0.2:1,'fontsize',26)
ylabel('contrast response scaling','fontsize',24)
xlabel('receptive field input strength','fontsize',24)

myFun=fittype('a*x');
rfFit=fit(scales,rfInputStrength',myFun);
totalScaleFactor=rfFit.a;



figure
set(gcf,'position',[1870         469         496         452])
plot(rfInputStrength,scales*totalScaleFactor,'sk','markersize',15,'linewidth',2)
axis([0 1.1 0 1.1])
line([0,1.1],[0,1.1],'linestyle','--','color','k')

title('cell ID 161, Fig6G, norm fit','fontsize',18)
set(gca,'xtick',0:0.2:1,'fontsize',26)
set(gca,'ytick',0:0.2:1,'fontsize',26)
ylabel('contrast response scaling','fontsize',24)
xlabel('receptive field input strength','fontsize',24)

%% cell 161 (3586), off midget - compilation

mySingleConeCell=my_cells(my_cells(:,2)==3586,1); % cell ID 3586, cell index 117 in datarun (single cones), 161 in datarunA (sta)
spikes=datarun.spikes{117};
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
myMaps=[11 13:19]; % starting from 0!
% myMaps=[11:19];
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
imagesc(-maps_comb(232:292,355:415))
title('Cell ID 3586')


% RF input strength

%plot sta of this cell
sta=squeeze(datarunA.stas.stas{161});

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
imagesc(resizedTC(232:292,355:415))
mySTA=resizedTC;



mySTA=imresize(sta(:,:,5),2, 'method','nearest');
figure
colormap gray
imagesc(mySTA(232:292,355:415))

maxx=max(max(mySTA(232:292,355:415)));
minn=min(min(mySTA(232:292,355:415)));

 % starting from 0!
myConeWeight=[];
myMap=zeros(61,61);tmp2=0;
figure
colormap gray
for myCone=1:length(myMaps)
    if myMaps(myCone)<10
        map=load(['/Volumes/Analysis/2011-12-13-2/stimuli/2011-12-13-2_f08_allcones/map-000',int2str(myMaps(myCone)),'.txt']);
    else
        map=load(['/Volumes/Analysis/2011-12-13-2/stimuli/2011-12-13-2_f08_allcones/map-00',int2str(myMaps(myCone)),'.txt']);
    end
    
%     figure
%     imagesc(map)
    tmp=mySTA.*map;
    tmp=tmp(232:292,355:415);
    subplot(3,3,myCone)
    tmp1=tmp;tmp1(1,1)=minn;tmp1(1,2)=maxx;
    imagesc(tmp1);    
    tmp2=tmp2+tmp;
    myConeWeight(myCone)=sum(tmp(:));
    tmp(tmp~=0)=sum(tmp(:));
    aa(myCone)=nnz(tmp);
    myMap=myMap+tmp;

end
tmp2(1,1)=minn;tmp2(1,2)=maxx;
subplot(3,3,9)
imagesc(tmp2)

figure
colormap gray
imagesc(myMap)

myConeWeight=myConeWeight./aa;

rfInputStrength=myConeWeight/min(myConeWeight);

origPicture=imread('/Users/alexth/Desktop/Peters_data_replica/OrigPicture.png');
origPicture=double(origPicture)/255;

figure
imagesc(origPicture)
x01=[190 651]; % 0 and 1
y01=[581 120]; % 0 and 1 in reversed order

coord=[651 122; 581 220; 499 246; 480 311; 464 338; 449 345; 373 320; 235 493];
PeterRFInputStrength=(coord(:,1)-x01(1))/diff(x01); % x coord (RF input strength)
PeterConeInputStrength=-(y01(1)-coord(:,2))/diff(y01); % y coord (con einput strength)

PeterRFInputStrength(end:-1:1)'-sort(rfInputStrength)'

figure
plot(PeterRFInputStrength(end:-1:1),myMeanResponses,'-*')

rfWeightsByPeter=[0.1045
    0.3872
    0.5522
    0.6014
    0.6308
    0.6746
    0.8501
    1.0000];

PeterRFInputStrength(end:-1:1)-rfWeightsByPeter


    0.1072
    0.3837
    0.4927
    0.5991
    0.6729    
    0.7427
    0.8383
    1.0000
    
svn=double(imread('/Users/alexth/Desktop/Peters_data_replica/Fig6B_fromPeter_SVN.png'));
figure
colormap gray
imagesc(svn(:,:,1))
colors=[0 2; 57 79; 112 97; 124 137; 133 152; 143 156; 184 142; 255 234]/255;
colors=1-colors*(1 - 0.1072);


origPicture=imread('/Users/alexth/Desktop/Peters_data_replica/Fig6C.png');
origPicture=double(origPicture)/255;
figure
imagesc(origPicture)
colors=[0.102 0.102; 0.322 0.396; 0.447 0.365; 0.502 0.561; 0.541 0.627; 0.588 0.643; 0.757 0.584; 1 0.933]';
colors=colors-min(colors(:));
colors=1-colors;

figure
hold on
plot(colors(1,:),colors(2,:),'dr','markersize',15,'linewidth',2)
plot(PeterRFInputStrength,PeterConeInputStrength,'sk','markersize',15,'linewidth',2)
h=legend('Fig6C (colors)','Fig6G (final)')
set(h,'Fontsize',18)
axis([0 1.1 0 1.1])
line([0,1.1],[0,1.1],'linestyle','--','color','k')

title('cell ID 161, Fig6G vs Fig6C','fontsize',18)
set(gca,'xtick',0:0.2:1,'fontsize',26)
set(gca,'ytick',0:0.2:1,'fontsize',26)
ylabel('contrast response scaling','fontsize',24)
xlabel('receptive field input strength','fontsize',24)

text(PeterRFInputStrength,PeterConeInputStrength+0.04,{'1','2','3','4','5','6','7','8'},'Fontsize',16,'color','k')

text(colors(1,:),colors(2,:)-0.04,{'1','2','3','4','5','6','7','8'},'Fontsize',16,'color','r')

figure
hold on
plot(PeterRFInputStrength,colors(1,:),'sk','markersize',15,'linewidth',2)
axis([0 1.1 0 1.1])
line([0,1],[0,1],'linestyle','--','color','k')


figure
hold on
plot(PeterConeInputStrength,colors(2,:),'sk','markersize',15,'linewidth',2)
axis([0 1.1 0 1.1])
line([0,1],[0,1],'linestyle','--','color','k')



figure
plot(myMeanResponses,'-*')
hold on
plot(rfInputStrength*6,'-*k')
plot(PeterRFInputStrength(end:-1:1)*6,'-*m')


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


for cnt=1
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

title({'contrast-response function, cell ID 161',['SCALES: ', num2str(scales',3)]},'fontsize',18)
set(gca,'xtick',0,'xticklabel','0','fontsize',26)
set(gca,'ytick',0:2:6,'yticklabel',{'0','2','4','6'},'fontsize',26)
ylabel('# of spikes','fontsize',28)
xlabel('contrast (arbitrary units)','fontsize',28)





figure
set(gcf,'position',[1870         469         496         452])
plot(rfInputStrength,scales,'sk','markersize',15,'linewidth',2)
axis([0 1.1 0 1.1])
line([0,1.1],[0,1.1],'linestyle','--','color','k')

title('cell ID 161, Fig6G','fontsize',18)
set(gca,'xtick',0:0.2:1,'fontsize',26)
set(gca,'ytick',0:0.2:1,'fontsize',26)
ylabel('contrast response scaling','fontsize',24)
xlabel('receptive field input strength','fontsize',24)



%% cell 61, off midget - compilation

mySingleConeCell=my_cells(my_cells(:,2)==1351,1); % cell ID 1906, cell index 48 in datarun
spikes=datarun.spikes{48};
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
    tt=map+map_tmp;
    imagesc(tt(210:280,100:165))
    title(int2str(i))

end
myMaps=[3:8 10:19]; % starting from 0!

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

figure
plot(myMeanResponses','-*')


% cone map with cone strength derived from relative response strength (NOT scaling)
myNormResponse=myMeanResponses/max(myMeanResponses);
maps_comb=0;
for i=1:length(myMaps)
    if myMaps(i)<10
        map=load(['/Volumes/Analysis/2011-12-13-2/stimuli/2011-12-13-2_f08_allcones/map-000',int2str(myMaps(i)),'.txt']);
    else
        map=load(['/Volumes/Analysis/2011-12-13-2/stimuli/2011-12-13-2_f08_allcones/map-00',int2str(myMaps(i)),'.txt']);
    end
    map(map>0)=myNormResponse(i);
    maps_comb=maps_comb+map;
end
figure
colormap gray
imagesc(-maps_comb(214:274,105:165))
title('strongest contrast response, not SCALING! Cell ID 61')


% RF input strength

%plot sta of this cell
sta=squeeze(datarunA.stas.stas{61});
mySTA=imresize(sta(:,:,5),2, 'method','nearest');
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
    tmp=tmp(214:274,105:165);
    myConeWeight(myCone)=sum(tmp(:));
    tmp(tmp~=0)=sum(tmp(:));
    myMap=myMap+tmp;

end

figure
colormap gray
imagesc(-myMap)
title('WN weights, cell ID 61')

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
[firstFit,gof]=fit(xdata,tmp,gausCum,'Start',[1.3,0.7,0.5],'Upper',[2, 1, 1]);


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
    [firstFit,gof]=fit(xdata,tmp,gausCum,'Start',[height mu sigma],'Upper',[2, 1, 1]);
    gof.rsquare    
end




figure
set(gcf,'position',[1621 525 526 481])
hold on
x=0.25:0.25:1;
markers='+o^sxd...........';
for i=1:size(myMeanResponses,1)
    plot(scales(i)*-x,myNormResponses(:,i)*max(myMeanResponses(:)),[markers(size(myMeanResponses,1)-i+1),'k'],'markersize',15,'linewidth',2)
end
axis tight
a=get(gca,'Xlim');
axis([a(1)*1.1 0 0 max(myMeanResponses(:))*1.1])
hold on
x=-1:0.05:0;
plot(x,cdf('norm',-x,firstFit.mu,firstFit.sigma)*firstFit.height*max(myMeanResponses(:)),'k');

title({'contrast-response function, cell ID 61',['SCALES: ', num2str(scales',3)]},'fontsize',18)
set(gca,'xtick',0,'xticklabel','0','fontsize',26)
set(gca,'ytick',0:2:6,'yticklabel',{'0','2','4','6'},'fontsize',26)
ylabel('# of spikes','fontsize',28)
xlabel('contrast (arbitrary units)','fontsize',28)





figure
set(gcf,'position',[1870         469         496         452])
plot(rfInputStrength,scales,'sk','markersize',15,'linewidth',2)
axis([0 1.1 0 1.1])
line([0,1.1],[0,1.1],'linestyle','--','color','k')

title('cell ID 61, Fig6G','fontsize',18)
set(gca,'xtick',0:0.2:1,'fontsize',26)
set(gca,'ytick',0:0.2:1,'fontsize',26)
ylabel('contrast response scaling','fontsize',24)
xlabel('receptive field input strength','fontsize',24)

[~,b]=sort(scales,'descend');

figure
for i=1:length(myMaps)
    plot(myMeanResponses(i,:),'-*','color',repmat(1-scales(i),1,3))
    hold on
    for j=1:4
        text(j+.03,myMeanResponses(i,j)+0.03,int2str(b(i)))
    end
end
axis([0.9 4.1 0.4 3.1])





fig6G=imread('/Users/alexth/Desktop/Peters_data_replica/Fig6G_cellF.jpeg');
figure
imagesc(fig6G)
coord=[505 216; 500 225; 461 301; 428 296; 422 302; 371 328; 280 323; 354 367; 366 384; 339 385; 282 361;...
    237 381; 212 409; 169 419; 192 456; 182 479; 181 504; 188 508; 151 506];
x01=[160 505];
y01=[537 187];

PeterRFinput=(coord(:,1)-x01(1))/diff(x01);
PeterFlashWeight=(y01(1)-coord(:,2))/(-diff(y01));

figure
plot(PeterRFinput',PeterFlashWeight,'vk','markersize',15,'linewidth',2)
axis([-0.1 1.1 -0.1 1.1])
line([-0.1,1.1],[-0.1,1.1],'linestyle','--','color','k')
title('Cell F from Fig6G')



a=sort(PeterRFinput);
b=sort(rfInputStrength');
[a(4:end) b]



fig6F=imread('/Users/alexth/Desktop/Peters_data_replica/Fig6F.png');
figure
imagesc(fig6F)

Peter6Frf=[26 26 39 65 86 141 164 185 199 207 211 225 229 238 238 254]/255;
Peter6FFlash=[65 26 89 70 26 150 170 187 189 162 191 219 238 255 232 229]/255;

Peter6Frf=1-(Peter6Frf-min(Peter6Frf));
Peter6FFlash=1-(Peter6FFlash-min(Peter6FFlash));


figure
plot(Peter6Frf,Peter6FFlash,'vk','markersize',15,'linewidth',2)
axis([-0.1 1.1 -0.1 1.1])
line([-0.1,1.1],[-0.1,1.1],'linestyle','--','color','k')
title('Cell F from Fig 6F')

