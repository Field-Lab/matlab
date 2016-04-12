function model = setup_cutting_metric(model,pixelSz)

stixelSz = pixelSz*3;
coarse_grid =  model.gridSzX/ (stixelSz);
coarse_dim = coarse_grid * coarse_grid;
coarse_pixels = eye(coarse_dim);

% go from high res to low res
A_coarse=sparse(model.gridSzX*model.gridSzY,coarse_grid*coarse_grid);
for ipix = 1:coarse_grid*coarse_grid
    A_coarse(:,ipix) = reshape(repelem(reshape(coarse_pixels(:,ipix),coarse_grid,coarse_grid),stixelSz,stixelSz),model.gridSzX*model.gridSzY,1);
end

% high res map of sub-units
su_highres = sparse(model.gridSzX*model.gridSzY,model.nSU);
iidx= 1:model.nCones;
for isu=1:model.nSU
    cones = iidx(model.cone_su_idx==isu);
    
    cMap =zeros(model.gridSzX,model.gridSzY);
    for icone =cones
        cMap = cMap + model.cone_data{icone}.coneMap;
    end
    su_highres(:,isu) = cMap(:);
end

% low resolution map of sub-units
su_lowres = A_coarse'*su_highres;

%% plot low res pixels on top of high res ones

mask2 = logical(ones(coarse_grid,coarse_grid));
sta_dim1 = size(mask2,1);
sta_dim2 = size(mask2,2);
indexedframe = reshape(1:sta_dim1*sta_dim2,[sta_dim1,sta_dim2]);
masked_frame = indexedframe(logical(mask2));

figure;
u_spatial_log = su_lowres;
for ifilt=1:model.nSU
    u_spatial = u_spatial_log(:,ifilt);
     u_spatial = reshape_vector(u_spatial,masked_frame,indexedframe);
     
subplot(ceil(sqrt(model.nSU)),floor((sqrt(model.nSU)))+1,ifilt);
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


%%
model.A_coarse = A_coarse;
model.su_highres = su_highres;
model.su_lowres = su_lowres;
end


