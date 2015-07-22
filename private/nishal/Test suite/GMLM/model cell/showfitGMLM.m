function showfitGMLM(fitGMLM2,text,mask)

nFilters = length(fitGMLM2.Linear.filter);

sta_dim1 = size(mask,1);
sta_dim2 = size(mask,2);
indexedframe = reshape(1:sta_dim1*sta_dim2,[sta_dim1,sta_dim2]);
masked_frame = indexedframe(logical(mask));

figure;
for ifilt=1:nFilters
subplot(2,2,ifilt)
u_spatial = reshape_vector(fitGMLM2.Linear.filter{ifilt}(1:length(masked_frame)),masked_frame,indexedframe);
imagesc((u_spatial));
%caxis([-0.3,0.3]);
colormap gray
colorbar
title(sprintf('%s : %d',text,ifilt));
end
