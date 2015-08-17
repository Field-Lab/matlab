
location_of_git_repo='/home/vision/Nishal/matlab/';
addpath(genpath([location_of_git_repo '/private/nora']));

% Java library

addpath(('/Volumes/Lab/Users/bhaishahster/GITs/matlab/private/nishal/create_act2'));
javaaddpath('/Volumes/Lab/Development/vision7/Vision.app/Contents/Resources/Java/Vision.jar');
addpath(genpath('/Volumes/Lab/Users/bhaishahster/GITs/matlab/private/nishal/Test suite/GLM/'));
addpath(genpath('/Volumes/Lab/Users/bhaishahster/GITs/matlab/private/nishal/Test suite/GLM2/'));
addpath(genpath('/Volumes/Lab/Users/bhaishahster/GITs/matlab/private/nishal'));
addpath(genpath('/Volumes/Lab/Users/bhaishahster/GITs/matlab/code'));
addpath(genpath('/Volumes/Lab/Users/bhaishahster/GITs/matlab/private/nishal/create_act_2/'));
%% Dataset details


WN_datafile = '2005-04-26-0/data009-dum/data009/data009';
WN_datafile_short=WN_datafile;
movie_xml = 'RGB-20-1-0.48-33333';
wrong_xml = 'RGB-20-1-0.48-11111';
stim_length=1800;% 
%% Load stimuli


datarun=load_data(WN_datafile)
datarun=load_params(datarun)
datarun=load_sta(datarun);
datarun=load_neurons(datarun);

for cellID = [922,809,258,1752,675,1918,1670,1999,1945]%datarun.cell_types{celltype_n}.cell_ids%input('Enter Cell ID');
%     
%   threshold_wrongSTA=sta_thr_using_wrong_movie(wrong_xml,datarun,stim_length,cellID);
%  cell_params.thres = threshold_wrongSTA;

extract_movie_response3;
inSU=0;
szz =[1,1;
    1,2;
    2,2;
    2,2;
    2,3;
    2,3;
    3,3;
    3,3;
    3,3;
    3,4;
    3,4;
    3,4;
    4,4;
    4,4;
    4,4;
    4,4;
    4,5];

 for nSU=[7,9,10,11:17,8,6,5,4,3,2,1]
inSU=inSU+1;     
    cellID

 %% gmlm
 mov_use = maskedMovdd;
 binnedResponsesbigd = spksGen;
filteredStimDim=size(mov_use,1)
 su_log=zeros(filteredStimDim,filteredStimDim);
 fitGMLM_log=cell(50,1);
 for ifit=1:10
     ifit
     %% fit
     close all
interval=1;
  [fitGMLM,f_val(ifit)] = fitGMLM_MEL_EM_bias(binnedResponsesbigd,mov_use,filteredStimDim,nSU,interval); 
   fitGMLM_log{ifit} = fitGMLM;
 
 %% plot
 h= figure('Color','w');

      mask = totalMaskAccept;
sta_dim1 = size(mask,1);
sta_dim2 = size(mask,2);
indexedframe = reshape(1:sta_dim1*sta_dim2,[sta_dim1,sta_dim2]);
masked_frame = indexedframe(logical(mask));

%u_spatial_log = zeros(40,40,nSU);

W=zeros(length(masked_frame),nSU);
for ifilt=1:nSU

u_spatial = reshape_vector(fitGMLM.Linear.filter{ifilt}(1:length(masked_frame)),masked_frame,indexedframe);
%u_spatial = reshape_vector(WN_uSq{ifilt},masked_frame,indexedframe);

subplot(szz(nSU,1),szz(nSU,2),ifilt)

imagesc(repelem(u_spatial(x_coord,y_coord),20,20));
colormap gray
%colorbar
title(sprintf('gmlm Filter: %d',ifilt));
axis square
u_spatial = u_spatial.*totalMaskAccept;
%u_spatial_log(:,:,ifilt) = u_spatial;
W(:,ifilt) = fitGMLM.Linear.filter{ifilt}(1:length(masked_frame)); 
set(gca,'xTick',[]);
set(gca,'yTick',[]);

end

if ~exist(sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2005_04_26_0/data009/%s/detailed/CellID_%d/exp_gmlm/SU_%d','cells',cellID,nSU),'dir');
    mkdir(sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2005_04_26_0/data009/%s/detailed/CellID_%d/exp_gmlm/SU_%d','cells',cellID,nSU));
end


hgexport(h,sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2005_04_26_0/data009/%s/detailed/CellID_%d/exp_gmlm/SU_%d/Cell%d_gmlm_su_%d_fit_%d.eps','cells',cellID,nSU,cellID,nSU,ifit));
close(h)
%% make occurance histogram
 k_est=W';
 su=double((k_est)==repmat(max((k_est),[],1),[nSU,1]));
 su_log = su_log + su'*su;
 end
 %% make near-neighbour co-occurance map
%  h=figure;
%  ilist= 1:length(masked_frame);
%  for ix=x_coord
%     for iy=y_coord
%     i1 =ilist(indexedframe(ix,iy)==masked_frame);
%     if(~isempty(i1))
%    datarun.cell_types{celltype_n}.nam
%     i2 = ilist(indexedframe(ix,iy+1)==masked_frame);
%     
%     if(~isempty(i2))
%         if(su_log(i1,i2)>0)
%     plot([ix,ix],[iy,iy+1],'LineWidth',su_log(i1,i2));
%     hold on;
%      text(ix,iy+0.5,sprintf('%d',su_log(i1,i2)));
%     hold on;
%         end
%     end
%     
%     i2 = ilist(indexedframe(ix,iy-1)==masked_frame);
%     if(~isempty(i2))
%           if(su_log(i1,i2)>0)
%     plot([ix,ix],[iy,iy-1],'LineWidth',su_log(i1,i2));
%     hold on;
%        text(ix,iy-0.5,sprintf('%d',su_log(i1,i2)));
%     hold on;
%           end
%     end
%     
%     i2 = ilist(indexedframe(ix-1,iy)==masked_frame);
%     if(~isempty(i2))
%            if(su_log(i1,i2)>0)
%     plot([ix-1,ix],[iy,iy],'LineWidth',su_log(i1,i2));
%     hold on;
%        text(ix-0.5,iy,sprintf('%d',su_log(i1,i2)));
%     hold on;
%            end
%     end
%     
%     
%     i2 = ilist(indexedframe(ix+1,iy)==masked_frame);
%     if(~isempty(i2))
%            if(su_log(i1,i2)>0)
%     plot([ix+1,ix],[iy,iy],'LineWidth',su_log(i1,i2));
%     hold on;
%        text(ix+0.5,iy,sprintf('%d',su_log(i1,i2)));
%     hold on;
%            end
%     end
%     
%     end
%     
%     end
%  end
%  
%  hgexport(h,sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2015_03_09_2/data038/Off parasol/detailed_subset/CellID_%d/gmlm/Cell%d_gmlm_su_%d.eps',cellID,cellID,nSU));
save(sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2005_04_26_0/data009/%s/detailed/CellID_%d/exp_gmlm/Cell%d_exp_gmlm_su_%d.mat','cells',cellID,cellID,nSU),'su_log','fitGMLM_log','x_coord','y_coord','totalMaskAccept','ttf','f_val');

 end
 end
 