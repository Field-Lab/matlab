%% load cones

piece = '2013-10-10-1';

% BW data
run = 'data000';
datarun = load_data(['/Volumes/Analysis/' piece '/streamed/' run '/' run]);
opt=struct('verbose',1,'load_params',1,'load_neurons',1,'load_obvius_sta_fits',false,'load_ei',false);
datarun=load_data(datarun,opt);
datarun = load_params(datarun,struct('verbose',1));  
datarun = set_polarities(datarun);
datarun = load_sta(datarun,'load_sta',[]);
datarun = load_cones(datarun,'300.00');
cones = datarun.cones.centers*9; % scale to match dll dimensions

tmp = detect(find_cone_data(datarun), @(cd) (regexp(cd, '300')));
tmp = [single_cone_path tmp,'/'];
load([tmp 'mydll.mat']);
dll = fliplr(mydll(:,:,1));

% RGB data
run = 'data002';
datarunA = load_data(['/Volumes/Analysis/' piece '/streamed/' run '/' run]);
opt = struct('verbose',1,'load_params',1,'load_neurons',1,'load_obvius_sta_fits',false,'load_ei',false);
datarunA=load_data(datarunA,opt);
datarunA = load_params(datarunA,struct('verbose',1));  
datarunA = set_polarities(datarunA);
datarunA = load_sta(datarunA,'load_sta',[]);
datarunA = load_cones(datarunA,'30.00');
conesA = datarunA.cones.centers*9; % scale to match dll dimensions

tmp = detect(find_cone_data(datarunA), @(cd) (regexp(cd, '30')));
tmp = [single_cone_path tmp,'/'];
load([tmp 'mydll.mat']);
dllA = fliplr(mydll(:,:,1));

clear tmp opt run piece mydll

% compare cone maps (check for movement), scaled and flipped (NOTE Y axis direction is reversed to match imagesc)
figure
hold on
plot(size(dll,1)-cones(:,1),cones(:,2),'x','color',[0.7 0.2 0.5])
plot(size(dllA,1)-conesA(:,1),conesA(:,2),'+','color',[0.2 0.7 0.5])
legend('BW','RGB')
set(gca, 'YDir','reverse')

% compare dlls (check for movement)
dll_combined=zeros([size(dll),3]);
dll_combined(:,:,1)=dll/max(dll(:));
dll_combined(:,:,2)=dllA/max(dllA(:));
dll_combined(dll_combined<0)=0;
figure
imagesc(dll_combined)
% add cone centers (already scaled, here flipped and shifted to match dll)
hold on
plot(size(dll,1)-cones(:,1),cones(:,2),'x','color',[0.7 0.2 0.5])
plot(size(dllA,1)-conesA(:,1),conesA(:,2),'+','color',[0.2 0.7 0.5])
legend('BW','RGB')




% load STA of an example cell in BW datarun
cellID=198;
datarun = load_sta(datarun,'load_sta',cellID);
sta=squeeze(datarun.stas.stas{datarun.cell_ids==cellID});
sta=sta(:,:,6);

% same cell in RGB datarun
cellID=197;
datarunA = load_sta(datarunA,'load_sta',cellID);
staA=datarunA.stas.stas{datarunA.cell_ids==cellID};
staA=staA(:,:,:,6);
staA=sum(staA,3);

sta_combined=zeros([size(sta),3]);
sta_combined(:,:,1)=sta/min(sta(:));
sta_combined(:,:,2)=staA/min(staA(:));
sta_combined(sta_combined<0)=0;
figure
imagesc(sta_combined)

% flip stas to match dll (which is also flipped to match camera image),
% resize 9 times and normalize to range from 0 to 1

staScaled=imresize(sta,9);
staScaled=staScaled/min(staScaled(:));
staScaled(staScaled<0)=0;
staScaled=fliplr(staScaled);

staScaledA=imresize(staA,9);
staScaledA=staScaledA/min(staScaledA(:));
staScaledA(staScaledA<0)=0;
staScaledA=fliplr(staScaledA);

% plot scaled STAs
sta_combined=zeros([size(staScaled),3]);
sta_combined(:,:,1)=staScaled;
sta_combined(:,:,2)=staScaledA;
figure
imagesc(sta_combined)

% plot together dll and sta for BW
dllovr=zeros(3000,3000,3);
dllovr(1:size(dll,1),10:size(dll,2)+10-1,2)=dll;
dllovr(5:size(staScaled,1)+5-1,1:size(staScaled,2),1)=staScaled;
figure
imagesc(dllovr);

% plot together dll and sta for RGB
dllovr=zeros(3000,3000,3);
dllovr(1:size(dllA,1),10:size(dllA,2)+10-1,2)=dllA;
dllovr(5:size(staScaledA,1)+5-1,1:size(staScaledA,2),1)=staScaledA;
figure
imagesc(dllovr);

%% load imaged cones

array_cones=double(imread('/Volumes/Analysis/2013-10-10-1/spot/Image22.tif'));
figure
colormap gray
imagesc(array_cones)

% resize to match dll size
myArray=imresize(array_cones,3.28);
figure;imagesc(myArray)

%% plot all together for BW
corCoefY=166; % to shift imaged array away from centered position in order to align selected cell STA best way. Smaller number moves array down relative to dll
corCoefX=35; % smaller number shifts array to the right

myMap=zeros(4000,5500,3);
myMap((corCoefY:size(myArray,1))+9-corCoefY,(((11+corCoefX):size(myArray,2))+9)-19-corCoefX,3)=myArray(corCoefY:end,(11+corCoefX):end)/max(myArray(:));
myMap((1:size(dll,1))+435,(1:size(dll,2))+1510+9,1)=dll/max(dll(:));
myMap((1:size(staScaled,1))+435+4,(1:size(staScaled,2))+1510,2)=staScaled;

figure
imagesc(myMap)
hold on

% put * for electrode 14 (where selected cell for sta is)
plot(2360+corCoefX, 1095+corCoefY,'*y')

%% plot all together for RGB

myMap=zeros(4000,5500,3);
myMap((corCoefY:size(myArray,1))+9-corCoefY,(((11+corCoefX):size(myArray,2))+9)-19-corCoefX,3)=myArray(corCoefY:end,(11+corCoefX):end)/max(myArray(:));
myMap((1:size(dllA,1))+435,(1:size(dllA,2))+1510+9,1)=dllA/max(dllA(:));
myMap((1:size(staScaledA,1))+435+4,(1:size(staScaledA,2))+1510,2)=staScaledA;

figure
imagesc(myMap)
hold on

% put * for electrode 14 (where selected cell for sta is)
plot(2360+corCoefX, 1095+corCoefY,'*y')

%% plot together array image (blue channel), BW dll (red channel), RGB dll (green channel)

corCoefY=166; % to shift imaged array away from centered position in order to align selected cell STA best way. Smaller number moves array down relative to dll
corCoefX=35; % smaller number shifts array to the right

myMap=zeros(4000,5500,3);
myMap((corCoefY:size(myArray,1))+9-corCoefY,(((11+corCoefX):size(myArray,2))+9)-19-corCoefX,2)=myArray(corCoefY:end,(11+corCoefX):end)/max(myArray(:));
myMap((1:size(dll,1))+435,(1:size(dll,2))+1510+9,3)=dll/max(dll(:));
myMap((1:size(dllA,1))+435,(1:size(dllA,2))+1510+9,1)=dllA/max(dllA(:));

figure
imagesc(myMap)
title('BLUE: array, RED: BW(data000), GREEN: RGB(data002)')
hold on

% put * for electrode 14 (where selected cell for sta is)
plot(2360+corCoefX, 1095+corCoefY,'*y')

% plot cone centers
plot(size(dll,1)-cones(:,1)+1510+9,cones(:,2)+435,'x','color',[0.7 0.2 0.5])
plot(size(dllA,1)-conesA(:,1)+1510+9,conesA(:,2)+435,'+','color',[0.2 0.7 0.5])

% axis([2000, 2350, 800,1500])
axis([2200, 2550, 1100,1500])


%% identify blue cones

figure
imagesc(myMap)
title('BLUE: array, RED: BW(data000), GREEN: RGB(data002)')
hold on
plot(size(dll,1)-cones(:,1)+1510+9,cones(:,2)+435,'x','color',[0.7 0.2 0.5])
plot(size(dllA,1)-conesA(:,1)+1510+9,conesA(:,2)+435,'+','color',[0.2 0.7 0.5])

plot(size(dllA,1)-conesA(tmp=='S',1)+1510+9,conesA(datarunA.cones.types=='S',2)+435,'*y')


%% local movement detection

nnd=pdist2(cones, conesA);
nnd_RGBtoBW=min(nnd);
nnd_BWtoRGB=min(nnd');

figure
hist(nnd_BWtoRGB,5:5:100)
axis([0 100 0 Inf])

figure
hold on
plot(cones(:,1),cones(:,2),'x','color',[0.7 0.2 0.5])
plot(conesA(:,1),conesA(:,2),'+','color',[0.2 0.7 0.5])
set(gcf,'Units','normalized')
set(gca,'Units','normalized')
ax = axis;
ap = get(gca,'Position');
hold off

[ac,ic]=sort(nnd_BWtoRGB);
isneighb=zeros(size(cones));

for i=1:find(ac>10,1)
    plot([cones(ic(i),1),conesA(ind,1)],[cones(ic(i),2),conesA(ind,2)],'.r')
    hold on
end
for i=find(ac>10,1):find(ac>45,1)
    [val,ind]=min(nnd(ic(i),:));
%     plot([cones(ic(i),1),conesA(ind,1)],[cones(ic(i),2),conesA(ind,2)])
    
    xo = [cones(ic(i),1),conesA(ind,1)];
    yo = [cones(ic(i),2),conesA(ind,2)];
    xp = (xo-ax(1))/(ax(2)-ax(1))*ap(3)+ap(1);
    yp = (yo-ax(3))/(ax(4)-ax(3))*ap(4)+ap(2);
    ah=annotation('arrow',xp,yp,'Color','r','HeadLength',3,'HeadWidth',5);
    drawnow
 
end

%% show directions on dll

figure
imagesc(myMap)
title('BLUE: array, RED: BW(data000), GREEN: RGB(data002)')
hold on
plot(size(dll,1)-cones(:,1)+1510+9,cones(:,2)+435,'x','color',[0.7 0.2 0.5])
plot(size(dllA,1)-conesA(:,1)+1510+9,conesA(:,2)+435,'+','color',[0.2 0.7 0.5])

for i=1:find(ac>30,1)
    [val,ind]=min(nnd(ic(i),:));
    plot(size(dll,1)-[cones(ic(i),1),conesA(ind,1)]+1510+9,[cones(ic(i),2),conesA(ind,2)]+435,'y')

end
