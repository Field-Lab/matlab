On = load('/Volumes/Analysis/nishal/recons_WN_whole_ON2.mat');
On2 = load('/Volumes/Analysis/nishal/recons_WN_whole_ON_nobias.mat');
Off = load('/Volumes/Analysis/nishal/recons_WN_whole_OFF.mat');
Off2 = load('/Volumes/Analysis/nishal/recons_WN_whole_OFF_nobias.mat');
OnOff  = load('/Volumes/Analysis/nishal/recons_WN_whole_PARASOL.mat');
OnOff2 = load('/Volumes/Analysis/nishal/recons_WN_whole_OnOFF_nobias.mat');


%%
recons_idx=[18000:21000]; % Doubt 

h=fspecial('gaussian',4,2);


for idx=recons_idx
    idx
close all
figure

subplot(4,3,1);
imagesc(On.mov_recons_full(:,:,idx));
colormap gray
axis image
caxis([-0.5,0.5]);
colorbar
%title('Spatial')

subplot(4,3,2);
imagesc(Off.mov_recons_full(:,:,idx));
colormap gray
caxis([-0.5,0.5]);
colorbar
axis image



subplot(4,3,3);
imagesc(OnOff.mov_recons_full(:,:,idx));
colormap gray
axis image
caxis([-0.5,0.5]);
colorbar

subplot(4,3,4);
imagesc(On.mov_pred_full(:,:,idx));
colormap gray
caxis([-0.5,0.5]);
colorbar
axis image
a=On.mov_recons_full(:,:,idx);
b=On.mov_pred_full(:,:,idx);
co=corr(a(:),b(:));
title(sprintf('Corr: %f',co));

subplot(4,3,5);
imagesc(Off.mov_pred_full(:,:,idx));
colormap gray
axis image
caxis([-0.5,0.5]);
colorbar
a=Off.mov_recons_full(:,:,idx);
b=Off.mov_pred_full(:,:,idx);
co=corr(a(:),b(:));
title(sprintf('Corr: %f',co));

subplot(4,3,6)
imagesc(OnOff.mov_pred_full(:,:,idx));
colormap gray
axis image
caxis([-0.5,0.5]);
colorbar
a=OnOff.mov_recons_full(:,:,idx);
b=OnOff.mov_pred_full(:,:,idx);
co=corr(a(:),b(:));
title(sprintf('Corr: %f',co));

subplot(4,3,7)
ppxx=imfilter(On.mov_recons_full(:,:,idx),h,'same').*(On.mov_recons_full(:,:,idx)~=0);
imagesc(ppxx)
colormap gray
axis image
caxis([-0.5,0.5]);
colorbar
a=ppxx(:);
b=On.mov_pred_full(:,:,idx);
co=corr(a(:),b(:));
title(sprintf('Corr: %f',co));

subplot(4,3,8)
ppxx=imfilter(Off.mov_recons_full(:,:,idx),h,'same').*(Off.mov_recons_full(:,:,idx)~=0);
imagesc(ppxx)
colormap gray
axis image
caxis([-0.5,0.5]);
colorbar
a=ppxx(:);
b=Off.mov_pred_full(:,:,idx);
co=corr(a(:),b(:));
title(sprintf('Corr: %f',co));

subplot(4,3,9)
ppxx=imfilter(OnOff.mov_recons_full(:,:,idx),h,'same').*(OnOff.mov_recons_full(:,:,idx)~=0);
imagesc(ppxx)
colormap gray
axis image
caxis([-0.5,0.5]);
colorbar
a=ppxx(:);
b=OnOff.mov_pred_full(:,:,idx);
co=corr(a(:),b(:));
title(sprintf('Corr: %f',co));

subplot(4,3,10);
imagesc(On2.mov_pred_full(:,:,idx));
colormap gray
caxis([-0.5,0.5]);
colorbar
axis image
a=imfilter(On.mov_recons_full(:,:,idx),h,'same').*(On.mov_recons_full(:,:,idx)~=0);
b=On2.mov_pred_full(:,:,idx);
co=corr(a(:),b(:));
title(sprintf('Corr: %f',co));

subplot(4,3,11);
imagesc(Off2.mov_pred_full(:,:,idx));
colormap gray
axis image
caxis([-0.5,0.5]);
colorbar
a=imfilter(Off.mov_recons_full(:,:,idx),h,'same').*(Off.mov_recons_full(:,:,idx)~=0);
b=Off2.mov_pred_full(:,:,idx);
co=corr(a(:),b(:));
title(sprintf('Corr: %f',co));


subplot(4,3,12)
imagesc(OnOff2.mov_pred_full(:,:,idx));
colormap gray
axis image
caxis([-0.5,0.5]);
colorbar
a=imfilter(OnOff.mov_recons_full(:,:,idx),h,'same').*(OnOff.mov_recons_full(:,:,idx)~=0);
b=OnOff2.mov_pred_full(:,:,idx);
co=corr(a(:),b(:));
title(sprintf('Corr: %f',co));



pause

end


%% Common mask calculation
Onmask=double(On.mov_pred_full(:,:,19000)~=0);

Offmask=double(Off.mov_pred_full(:,:,19000)~=0);
OnOffCommonMask=Onmask.*Offmask;
OnOffmask=double(OnOff.mov_pred_full(:,:,19000)~=0);
[row,col]=find(OnOffCommonMask==1);
%%
numPix=1;

commPixIdx= randsample(length(row),numPix);
icnt=0;
ppxx=zeros(size(On.mov_recons_full));

h=fspecial('gaussian',4,2);
for idx=recons_idx
ppxx(:,:,idx)=imfilter(On.mov_recons_full(:,:,idx),h,'same').*(On.mov_recons_full(:,:,idx)~=0);
end

plot_idx=idx;

figure;
for iPixIdx=commPixIdx'
    
pixX=15%row(iPixIdx);
pixY=15%col(iPixIdx);
recons_idx=[18000:21000]; % Doubt 

icnt=icnt+1;
subplot(numPix,3,icnt)
a=squeeze(ppxx(pixX,pixY,recons_idx));
b=squeeze(On.mov_pred_full(pixX,pixY,recons_idx));
stairs(a,'r');
hold on
stairs(b,'b');
title(sprintf('Corr: %f , Cells %d',corr(a,b),On.filter_collection(pixX,pixY).sig_cells));
ylim([-1,1]);

icnt=icnt+1;
subplot(numPix,3,icnt)
a=squeeze(ppxx(pixX,pixY,recons_idx));
b=squeeze(Off.mov_pred_full(pixX,pixY,recons_idx));
stairs(a,'r');
hold on
stairs(b,'g');
title(sprintf('Corr: %f , Cells %d',corr(a,b),Off.filter_collection(pixX,pixY).sig_cells));
ylim([-1,1])

icnt=icnt+1;
subplot(numPix,3,icnt)
a=squeeze(ppxx(pixX,pixY,recons_idx));
b=squeeze(OnOff.mov_pred_full(pixX,pixY,recons_idx));
stairs(a,'r');
hold on
stairs(b,'k');
title(sprintf('Corr: %f , Cells %d',corr(a,b),OnOff.filter_collection(pixX,pixY).sig_cells));
ylim([-1,1])

end


%%
dataCorrCell=struct([]);
dataCorrCell(1).type='On';
dataCorrCell(2).type='Off';
dataCorrCell(3).type='OnOff';

for itype=1:3
dataCorrCell(itype).Corr=[];
dataCorrCell(itype).Cell=[];
end

for iPixIdx=1:length(row)
    iPixIdx
pixX=row(iPixIdx);
pixY=col(iPixIdx);
recons_idx=[18000:21000]; % Doubt 


a=squeeze(ppxx(pixX,pixY,recons_idx));
b=squeeze(On.mov_pred_full(pixX,pixY,recons_idx));
sprintf('On: Corr: %f , Cells %d',corr(a,b),On.filter_collection(pixX,pixY).sig_cells);
dataCorrCell(1).Corr(iPixIdx)=corr(a,b);
dataCorrCell(1).Cell(iPixIdx)=On.filter_collection(pixX,pixY).sig_cells;

a=squeeze(ppxx(pixX,pixY,recons_idx));
b=squeeze(Off.mov_pred_full(pixX,pixY,recons_idx));
sprintf('Off: Corr: %f , Cells %d',corr(a,b),Off.filter_collection(pixX,pixY).sig_cells);
dataCorrCell(2).Corr(iPixIdx)=corr(a,b);
dataCorrCell(2).Cell(iPixIdx)=Off.filter_collection(pixX,pixY).sig_cells;


a=squeeze(ppxx(pixX,pixY,recons_idx));
b=squeeze(OnOff.mov_pred_full(pixX,pixY,recons_idx));
sprintf('On Off: Corr: %f , Cells %d',corr(a,b),Off.filter_collection(pixX,pixY).sig_cells);
dataCorrCell(3).Corr(iPixIdx)=corr(a,b);
dataCorrCell(3).Cell(iPixIdx)=OnOff.filter_collection(pixX,pixY).sig_cells;

end

figure;
scatter(dataCorrCell(1).Cell,dataCorrCell(1).Corr,40,'r');
hold on
scatter(dataCorrCell(2).Cell,dataCorrCell(2).Corr,30,'g');
hold on
scatter(dataCorrCell(3).Cell,dataCorrCell(3).Corr,20,'b');
legend('On','Off','On Off');
ylabel('Correlation');
xlabel('Number of Cells');
hold on 
line(1:10,0*[1:10]);

figure;
subplot(1,2,1);
scatter(dataCorrCell(1).Cell,dataCorrCell(3).Corr);
subplot(1,2,2);
scatter(dataCorrCell(2).Cell,dataCorrCell(3).Corr);
% 
% figure;
% surface(dataCorrCell(1).Cell,dataCorrCell(2).Cell,dataCorrCell(3).Corr);


CorrMat=zeros(max(dataCorrCell(1).Cell),max(dataCorrCell(2).Cell));
for iOnCells=1:size(CorrMat,1)
    for iOffCells=1:size(CorrMat,2)
    indices=(dataCorrCell(1).Cell==iOnCells)&(dataCorrCell(2).Cell==iOffCells);
    CorrMat(iOnCells,iOffCells)=mean(dataCorrCell(3).Corr(indices));
    end
end
figure;
surf(CorrMat)
xlabel('On Cells')
ylabel('Off Cells')
zlabel('Mean Correlation')




