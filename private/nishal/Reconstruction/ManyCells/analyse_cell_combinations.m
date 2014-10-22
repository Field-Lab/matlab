On = load('/Volumes/Analysis/nishal/recons_WN_whole_ON2.mat');
Off = load('/Volumes/Analysis/nishal/recons_WN_whole_OFF.mat');
OnOff  = load('/Volumes/Analysis/nishal/recons_WN_whole_PARASOL.mat');
recons_idx=[18000:21000]; % Doubt 

h=fspecial('gaussian',4,2);

figure
for idx=recons_idx
    idx
subplot(3,3,1);
imagesc(On.mov_recons_full(:,:,idx));
colormap gray
axis image
caxis([-0.5,0.5]);
colorbar
%title('Spatial')

subplot(3,3,2);
imagesc(Off.mov_recons_full(:,:,idx));
colormap gray
caxis([-0.5,0.5]);
colorbar
axis image



subplot(3,3,3);
imagesc(OnOff.mov_recons_full(:,:,idx));
colormap gray
axis image
caxis([-0.5,0.5]);
colorbar

subplot(3,3,4);
imagesc(On.mov_pred_full(:,:,idx));
colormap gray
caxis([-0.5,0.5]);
colorbar
axis image
a=On.mov_recons_full(:,:,idx);
b=On.mov_pred_full(:,:,idx);
co=corr(a(:),b(:));
title(sprintf('Corr: %f',co));

subplot(3,3,5);
imagesc(Off.mov_pred_full(:,:,idx));
colormap gray
axis image
caxis([-0.5,0.5]);
colorbar
a=Off.mov_recons_full(:,:,idx);
b=Off.mov_pred_full(:,:,idx);
co=corr(a(:),b(:));
title(sprintf('Corr: %f',co));

subplot(3,3,6)
imagesc(OnOff.mov_pred_full(:,:,idx));
colormap gray
axis image
caxis([-0.5,0.5]);
colorbar
a=OnOff.mov_recons_full(:,:,idx);
b=OnOff.mov_pred_full(:,:,idx);
co=corr(a(:),b(:));
title(sprintf('Corr: %f',co));

subplot(3,3,7)
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

subplot(3,3,8)
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

subplot(3,3,9)
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


pause

end


