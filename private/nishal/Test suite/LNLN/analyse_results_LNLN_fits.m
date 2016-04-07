fits = load('/Volumes/Lab/Users/bhaishahster/GITs/python_env/modelLNLN_result_small.mat');
mdl = load('/Volumes/Lab/Users/bhaishahster/GITs/python_env/modelLNLN_small2.mat');
model=mdl.model;
mask2 = zeros(10,10);
mask2(1:5,1:5)=1;
mask2=logical(mask2);

% iSU=4;
% idx = 1+(iSU-1)*3;
% K = fits.op_list{idx};

iSU=15
K = fits.K_hierarchy
nSU= iSU


sta_dim1 = size(mask2,1);
sta_dim2 = size(mask2,2);
indexedframe = reshape(1:sta_dim1*sta_dim2,[sta_dim1,sta_dim2]);
masked_frame = indexedframe(logical(mask2));

nSU=iSU
for ifilt=1:nSU
%subplot(szz(nSU,1),szz(nSU,2),ifilt)
subplot(1,nSU,ifilt);
%u_spatial = reshape_vector(fitGMLM{1}.Linear.filter{ifilt}(1:length(masked_frame)),masked_frame,indexedframe);
%u_spatial = -reshape(K(:,ifilt),10,10);
 u_spatial = -reshape_vector(K(:,ifilt),masked_frame,indexedframe);
strongestFrame = u_spatial/max(abs(u_spatial(:)))+0.5;
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

%print(h,'-dpdf',sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/modelCellLNLN/model1/Exp_ASM_MEL_EM_filters SU%d.pdf',nSU));

% plot weights
figure;
plot(model.SU_gang_weights.*sqrt(sum(model.cone_to_SU_connection,2)),'*')
title('SU to ganglion Weight distribution');


%%
su_strength = model.SU_gang_weights.*sqrt(sum(model.cone_to_SU_connection,2));
true_SU_legend=cell(4,1);
for itruesu = 1:4
true_SU_legend{itruesu} = sprintf('true su %d, mag: %0.02f',itruesu,su_strength(itruesu));
end

f=figure
for iSU=4;
    idx = 1+(iSU-1)*3;
    %K = fits.op_list{idx};
    %K=fits.K_ASM6;
    K = fits.K_hierarchy
    ex_SU_legend=cell(iSU,1);
    for iexsu = 1:iSU
         ex_SU_legend{iexsu} = sprintf('extracted su %d',iexsu);
    end
    % radar chart
    nSU = iSU
    dot_filter_su=zeros(nSU,model.nSU);
    for ifilter = 1:size(K,2)
        %     u_spatial = reshape_vector(fitGMLM{1}.Linear.filter{ifilter}(1:length(masked_frame)),masked_frame,indexedframe);
         u_spatial = reshape_vector(K(:,ifilter),masked_frame,indexedframe);
        %u_spatial = reshape(K(:,ifilter),10,10);
        strongestFrame = u_spatial;
        szstr = size(strongestFrame,1);
        ssf = repelem(strongestFrame,model.gridSzX/szstr,model.gridSzX/szstr,1);
        [SU_inp]=filter_su_dot_LNLN(model,ssf);
        
        dot_filter_su(ifilter,:) = SU_inp'./sqrt(sum(model.cone_to_SU_connection.^2,2))'; % divide by sqrt(#cones) ?
    end
    h =spider(dot_filter_su(4:7,:)',sprintf('Extracted filters: %d',nSU),[],true_SU_legend,ex_SU_legend,f.Number);
    
end