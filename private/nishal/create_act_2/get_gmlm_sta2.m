
function [stas_true,stas_rand] = get_gmlm_sta2(WN_datafile,movie_xml,stim_length,cellID_list,nSU_list,destination,sta_depth,contrast_factor)


%% Load stimuli


datarun=load_data(WN_datafile)
datarun=load_params(datarun)
datarun=load_sta(datarun);
datarun=load_neurons(datarun);

icell=0;

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

for nSU = nSU_list%1:2 %1:10%1:filteredStimDim
    nSU
 binnedResponsesbigd = spksGen(trainData);
binnedResponsesbigd_hr = spksGen_hr(trainData_hr);
 mov_use=maskedMovdd(:,trainData);
  filteredStimDim=size(mov_use,1); 
 
 [fitGMLM,output] = fitGMLM_MEL_EM(binnedResponsesbigd,mov_use,filteredStimDim,nSU,interval);  
 %[fitGMLM,f_val(nSU)] = fitGMLM_MEL_EM_bias(binnedResponsesbigd,mov_use,filteredStimDim,nSU,interval); 
% [fitGMLM,output] = fitGMLM_MEL_EM_power2(binnedResponsesbigd,mov_use,filteredStimDim,nSU,interval,2);  
 % fitGMLM_log(ifit).fitGMLM = fitGMLM;  %%
   [fitGMLM,output]= fitGMLM_full(fitGMLM,binnedResponsesbigd_hr,mov_use);
%   fitGMLM_full2_log{nSU}=fitGMLM;
%   figure;
%   plot(fitGMLM.hist.hexpanded)
fitGMLM_log{nSU} = fitGMLM;

close all;
end

 
 %save(sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2015_03_09_2/data031/Cell%d_full',cellID),'fitGMLM_log','fitGMLM_full2_log','mov_use','binnedResponsesbigd','filteredStimDim','interval','totalMaskAccept2','totalMaskAccept','x_coord','y_coord');
%% Compute STC 
%[WNSTA,WNSTC,WN_uSq]=compute_STA_STC(binnedResponsesbigd,mov_use);


 
  %% Show learned filters;


mask = totalMaskAccept;
sta_dim1 = size(mask,1);
sta_dim2 = size(mask,2);
indexedframe = reshape(1:sta_dim1*sta_dim2,[sta_dim1,sta_dim2]);
masked_frame = indexedframe(logical(mask));

for nSU = nSU_list
    close all;
    figure;
    fitGMLM=fitGMLM_log{nSU};
u_spatial_log = zeros(sta_dim1,sta_dim2,nSU);

for ifilt=1:nSU

u_spatial = reshape_vector(fitGMLM.Linear.filter{ifilt}(1:length(masked_frame)),masked_frame,indexedframe);
%u_spatial = reshape_vector(WN_uSq{ifilt},masked_frame,indexedframe);

subplot(ceil(floor(sqrt(nSU))),ceil(nSU/(floor(sqrt(nSU)))),ifilt);
imagesc(u_spatial(x_coord,y_coord));
colormap gray
%colorbar
%title(sprintf('ki: %d, e(b): %0.02f',ifilt,exp(fitGMLM.Linear.bias{ifilt})),'FontSize',10);
axis image
set(gca,'xTick',[]);set(gca,'yTick',[]);
u_spatial = u_spatial.*totalMaskAccept;
u_spatial_log(:,:,ifilt) = u_spatial;
end
suptitle(sprintf('nSU: %d',nSU));

sus_elim = [];%input('Sub-unit indices which are bad (1-..)?');good_su_idx=ones(nSU,1);
good_su_idx = ones(nSU,1);
if(~isempty(sus_elim))
good_su_idx(sus_elim)=0;
end
fitGMLM_log{nSU}.good_su_idx=good_su_idx;
end

pause(0.3);
%% make STAs
ttf=squeeze(ttf);
for nSU=nSU_list
ista=1; jsta=1;stas_true=cell(nSU,1);stas_rand=cell(nSU,1);
for isu=1:nSU
    u_st = repmat(u_spatial_log(:,:,isu)',[1,1,3,30]);
    for itime=1:30
        u_st(:,:,:,end-itime+1) = u_st(:,:,:,end-itime+1)*double(ttf(itime));
    end
    stas_true{ista}=u_st;
    ista=ista+1;
end

fitGMLM_log{nSU}.stas_true=stas_true;
urand =zeros(sta_dim1,sta_dim2,nSU);
for ix=1:sta_dim1
    for iy=1:sta_dim2
    p=randperm(nSU);
    urand(ix,iy,:)=u_spatial_log(ix,iy,p);
    end
end

for jsu=1:nSU
    u_st = repmat(urand(:,:,jsu)',[1,1,3,30]);
    for itime=1:30
        u_st(:,:,:,end-itime+1) = u_st(:,:,:,end-itime+1)*ttf(itime);
    end
    stas_rand{jsta}=u_st;
    jsta=jsta+1;
end

fitGMLM_log{nSU}.stas_rand=stas_rand;
end
%% save
if(~isdir(sprintf('/Volumes/Lab/Users/bhaishahster/%s',destination)))
    mkdir(sprintf('/Volumes/Lab/Users/bhaishahster/%s',destination));
end

save(sprintf('/Volumes/Lab/Users/bhaishahster/%s/Cell_%d.mat',destination,cellID),'fitGMLM_log','cellID','mask','x_coord','y_coord','ttf','-v7.3');
end


end
