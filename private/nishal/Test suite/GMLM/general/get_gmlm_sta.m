
function [stas_true,stas_rand] = get_gmlm_sta(WN_datafile,movie_xml,stim_length,cellID_list,nSU_list,destination,sta_depth)


%% Load stimuli


datarun=load_data(WN_datafile)
datarun=load_params(datarun)
datarun=load_sta(datarun);
datarun=load_neurons(datarun);

icell=0;
ista=1; jsta=1;
user_STA_depth = sta_depth;
for cellID = cellID_list %[datarun.cell_types{12}.cell_ids,datarun.cell_types{1}.cell_ids,datarun.cell_types{2}.cell_ids];
    cellID
icell=icell+1;
    %cell_glm_fit = sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_03_09_2/data038/CellType_OFF parasol/CellID_%d.mat',cellID);
%load(cell_glm_fit);
%%

extract_movie_response2;

 %% EM like Max Expected Likelihood .. 
 interval=1;
%  idx = [1:end]
trainData=[1:floor(length(spksGen))];
trainData_hr = [1:floor(length(spksGen_hr))];


 nSU = nSU_list(icell)%1:2 %1:10%1:filteredStimDim
 binnedResponsesbigd = spksGen(trainData);
binnedResponsesbigd_hr = spksGen_hr(trainData_hr);
 mov_use=maskedMovdd(:,trainData);
  filteredStimDim=size(mov_use,1); 
 
% [fitGMLM,output] = fitGMLM_MEL_EM(binnedResponsesbigd,mov_use,filteredStimDim,nSU,interval);  
 [fitGMLM,f_val(nSU)] = fitGMLM_MEL_EM_bias(binnedResponsesbigd,mov_use,filteredStimDim,nSU,interval); 
 fitGMLM_log{nSU} = fitGMLM;
 %[fitGMLM,output] = fitGMLM_MEL_EM_power2(binnedResponsesbigd,mov_use,filteredStimDim,nSU,interval,2);  
 % fitGMLM_log(ifit).fitGMLM = fitGMLM;  %%
%   [fitGMLM,output]= fitGMLM_full(fitGMLM,binnedResponsesbigd_hr,mov_use);
%   fitGMLM_full2_log{nSU}=fitGMLM;
%   figure;
%   plot(fitGMLM.hist.hexpanded)
 

 
 %save(sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2015_03_09_2/data031/Cell%d_full',cellID),'fitGMLM_log','fitGMLM_full2_log','mov_use','binnedResponsesbigd','filteredStimDim','interval','totalMaskAccept2','totalMaskAccept','x_coord','y_coord');
%% Compute STC 
%[WNSTA,WNSTC,WN_uSq]=compute_STA_STC(binnedResponsesbigd,mov_use);


 
  %% Show learned filters;
figure;

    
 %   fitGMLM_full2 = fitGMLM_full2_log{nSU};

mask = totalMaskAccept;
sta_dim1 = size(mask,1);
sta_dim2 = size(mask,2);
indexedframe = reshape(1:sta_dim1*sta_dim2,[sta_dim1,sta_dim2]);
masked_frame = indexedframe(logical(mask));

u_spatial_log = zeros(sta_dim1,sta_dim2,nSU);

for ifilt=1:nSU

u_spatial = reshape_vector(fitGMLM.Linear.filter{ifilt}(1:length(masked_frame)),masked_frame,indexedframe);
%u_spatial = reshape_vector(WN_uSq{ifilt},masked_frame,indexedframe);

subplot(7,2,ifilt)
imagesc(u_spatial(x_coord,y_coord));
colormap gray
colorbar
title(sprintf('GMLM Filter: %d',ifilt));
axis image
u_spatial = u_spatial.*totalMaskAccept;
u_spatial_log(:,:,ifilt) = u_spatial;
end

% [v,I] = max(u_spatial_log,[],3);
% 
% xx=I.*(v>0.2);
% xx=xx(x_coord,y_coord);
% figure;
% % subplot(2,2,nSU);
%  imagesc(xx);
%  title(sprintf('Num SU : %d',nSU));
pause(0.3);
%% make STAs
ttf=squeeze(ttf);
% figure;
for isu=1:nSU
    u_st = repmat(u_spatial_log(:,:,isu)',[1,1,3,30]);
    for itime=1:30
        u_st(:,:,:,end-itime+1) = u_st(:,:,:,end-itime+1)*double(ttf(itime));
    end
    stas_true{ista}=u_st;
    
%     for itime=1:30
%     imagesc(stas_true{ista}(:,:,itime));
%     caxis([min(stas_true{ista}(:)),max(stas_true{ista}(:))]);
%     colorbar
%     colormap gray
%     axis image
%     pause(1/120)
%     end
    ista=ista+1;
    
  
end



urand =zeros(sta_dim1,sta_dim2,nSU);
for ix=1:sta_dim1
    for iy=1:sta_dim2
    p=randperm(nSU);
    urand(ix,iy,:)=u_spatial_log(ix,iy,p);
    end
end
% 
% figure;

for jsu=1:nSU
    u_st = repmat(urand(:,:,jsu)',[1,1,3,30]);
    for itime=1:30
        u_st(:,:,:,end-itime+1) = u_st(:,:,:,end-itime+1)*ttf(itime);
    end
    stas_rand{jsta}=u_st;
%     
%     for itime=1:30
%     imagesc(stas_rand{jsta}(:,:,itime));
%     caxis([min(stas_rand{jsta}(:)),max(stas_rand{jsta}(:))]);
%     colorbar
%     colormap gray
%     axis image
%     pause(1/120)
%     end
    jsta=jsta+1;
end

%% save
if(~isdir(sprintf('/Volumes/Lab/Users/bhaishahster/%s',destination)))
    mkdir(sprintf('/Volumes/Lab/Users/bhaishahster/%s',destination));
end

save(sprintf('/Volumes/Lab/Users/bhaishahster/%s/Cell_%d.mat',destination,cellID),'fitGMLM','cellID','stas_true','stas_rand','mask','x_coord','y_coord','ttf');
end


end
