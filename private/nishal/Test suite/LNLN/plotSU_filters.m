function plotSU_filters (u_spatial_log,model)
figure;
mask2=logical(ones(10,10));
sta_dim1 = size(mask2,1);
sta_dim2 = size(mask2,2);
indexedframe = reshape(1:sta_dim1*sta_dim2,[sta_dim1,sta_dim2]);
masked_frame = indexedframe(logical(mask2));
nSU = size(u_spatial_log,2);

toplot = u_spatial_log;
for ifilt=1:nSU
    u_spatial = toplot(:,ifilt);
     u_spatial = reshape_vector(u_spatial,masked_frame,indexedframe);
subplot(1,nSU,ifilt);
strongestFrame = -1*u_spatial/max(abs(u_spatial(:)))+0.5;
szstr = size(strongestFrame,1);
ssf = repelem(strongestFrame,model.gridSzX/szstr,model.gridSzX/szstr,3);
mmask = sum(model.totalConeMap3D,3)==0;
xx = repelem(mmask,1,1,3).*ssf; 
aa = model.totalConeMap3D;  aa(repelem(sum(aa,3),1,1,3)<0.5)=aa(repelem(sum(aa,3)<0.5,1,1,3))+0.5;aa =(aa-0.5)/max(aa(:)-0.5) + 0.5;
mxt = ((aa(repelem(sum(aa,3)~=1.5,1,1,3))));mx= max(mxt(:));
xx(repelem(sum(aa,3)~=1.5,1,1,3))= mxt;
xx = xx(1:ceil(max(model.conesX))+40,1:ceil(max(model.conesY))+40,:);

% mag = repelem(sum(model.totalConeMap3D,3)==0,1,1,3).*ssf + model.totalConeMap3D*2;
% mag = mag(1:ceil(max(model.conesX))+40,1:ceil(max(model.conesY))+40,:);

imagesc(xx);
axis image
set(gca,'xTick',[]);
set(gca,'yTick',[]);

title(sprintf('SU # %d',ifilt));
%caxis([-max(mag(:)),max(mag(:))]);
%caxis([0,max(aa(:))])
end



end