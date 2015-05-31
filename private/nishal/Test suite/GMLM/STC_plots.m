%% GMLM with MEL and EM for normal cells at coarser resolution .. 

location_of_git_repo='/home/vision/Nishal/matlab/';
addpath(genpath([location_of_git_repo '/private/nora']));

% Java library

addpath(('~/Nishal/matlab/private/nishal/create_act2'));
javaaddpath('/Volumes/Lab/Development/vision7/Vision.app/Contents/Resources/Java/Vision.jar');
addpath(genpath('~/Nishal/matlab/private/nishal/Test suite/GLM/'));
addpath(genpath('~/Nishal/matlab/private/nishal/Test suite/GLM2/'));
addpath(genpath('~/Nishal/matlab/private/nishal'));
addpath(genpath('~/Nishal/matlab/code'));
addpath(genpath('~/Nishal/matlab/private/nishal/create_act_2/'));
%% Dataset details
% WN_datafile = '2015-03-09-2/data031/data031';
% WN_datafile_short='2015-03-09-2/data031/data031';
% movie_xml = 'BW-8-6-0.48-11111-40x40';
% stim_length=1800;% in seconds

WN_datafile = '2015-03-09-2/streamed/data038/data038';
WN_datafile_short='2015-03-09-2/streamed/data038/data038';
movie_xml = 'BW-8-2-0.48-11111-40x40';
stim_length=1800;% 
%% Load stimuli


datarun=load_data(WN_datafile)
datarun=load_params(datarun)
datarun=load_sta(datarun);
datarun=load_neurons(datarun);

for cellID = datarun.cell_types{2}.cell_ids %[datarun.cell_types{12}.cell_ids,datarun.cell_types{1}.cell_ids,datarun.cell_types{2}.cell_ids];
    cellID
%cell_glm_fit = sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_03_09_2/data038/CellType_OFF parasol/CellID_%d.mat',cellID);
%load(cell_glm_fit);
%%
extract_movie_response2;

 %% EM like Max Expected Likelihood .. 
 interval=1;
%  idx = [1:end]
trainData=[1:floor(length(spksGen)*0.8)];
trainData_hr = [1:floor(length(spksGen_hr)*0.8)];

testData=[floor(length(spksGen)*0.8)+1 :length(spksGen)];
testData_hr=[floor(length(spksGen_hr)*0.8)+1 :length(spksGen_hr)];

 binnedResponsesbigd = spksGen(trainData);
binnedResponsesbigd_hr = spksGen_hr(trainData_hr);
 mov_use=maskedMovdd(:,trainData);
  filteredStimDim=size(mov_use,1); 
 
 
 
 %save(sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2015_03_09_2/data031/Cell%d_full',cellID),'fitGMLM_log','fitGMLM_full2_log','mov_use','binnedResponsesbigd','filteredStimDim','interval','totalMaskAccept2','totalMaskAccept','x_coord','y_coord');
%% Compute STC 
[WNSTA,WNSTC,WN_uSq,h]=compute_STA_STC(binnedResponsesbigd,mov_use);

  if(~isdir(sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2015_03_09_2/data038/Off parasol/CellID_%d',cellID)))
        mkdir(sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2015_03_09_2/data038/Off parasol/CellID_%d',cellID));
  end
  hgexport(h,sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2015_03_09_2/data038/Off parasol/CellID_%d/spectrum.eps',cellID));



 
%% Prediction 
 
  %% Show learned filters;
h=figure;
 nSU=3
    
 %   fitGMLM_full2 = fitGMLM_full2_log{nSU};
      mask = totalMaskAccept;
sta_dim1 = size(mask,1);
sta_dim2 = size(mask,2);
indexedframe = reshape(1:sta_dim1*sta_dim2,[sta_dim1,sta_dim2]);
masked_frame = indexedframe(logical(mask));

u_spatial_log = zeros(40,40,nSU);
u_spatial = reshape_vector(WN_uSq{1},masked_frame,indexedframe);

subplot(2,2,1)
imagesc(repelem(u_spatial(x_coord,y_coord),10,10));
colormap gray
colorbar
title(('STA'));
axis square
u_spatial = u_spatial.*totalMaskAccept;
u_spatial_log(:,:,ifilt) = u_spatial;


for ifilt=1:nSU

%u_spatial = reshape_vector(fitGMLM_full2.Linear.filter{ifilt}(1:length(masked_frame)),masked_frame,indexedframe);
u_spatial = reshape_vector(WN_uSq{ifilt+1},masked_frame,indexedframe);

subplot(2,2,ifilt+1)
imagesc(repelem(u_spatial(x_coord,y_coord),10,10));
colormap gray
colorbar
title(sprintf('STC Filter: %d',ifilt));
axis square
u_spatial = u_spatial.*totalMaskAccept;
u_spatial_log(:,:,ifilt) = u_spatial;
end

[v,I] = max(u_spatial_log,[],3);

xx=I.*(v>0.2);
xx=xx(x_coord,y_coord);
%figure;
% subplot(2,2,nSU);
% imagesc(xx);
% title(sprintf('Num SU : %d',nSU));

  if(~isdir(sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2015_03_09_2/data038/Off parasol/CellID_%d',cellID)))
        mkdir(sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2015_03_09_2/data038/Off parasol/CellID_%d',cellID));
  end
  hgexport(h,sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2015_03_09_2/data038/Off parasol/CellID_%d/STA_STC%d.eps',cellID,nSU));


end
