
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

WN_datafile = '2008-08-27-6/data009/data009';
WN_datafile_short=WN_datafile;
movie_xml = 'RGB-10-1-0.48-11111';
stim_length=1800;% 
%% Load stimuli


datarun=load_data(WN_datafile)
datarun=load_params(datarun)
datarun=load_sta(datarun);
datarun=load_neurons(datarun);

icell_l=0;
cells = [datarun.cell_types{1}.cell_ids,datarun.cell_types{8}.cell_ids];
ttf_log = zeros(30,length(cells));
totalMaskAccept_log = zeros(64*32,length(cells));
Y = zeros(length(cells),216000);

for cellID = cells
    icell_l=icell_l+1
extract_movie_response4;
ttf_log(:,icell_l) =ttf;
totalMaskAccept_log(:,icell_l) = totalMaskAccept(:);
Y(icell_l,:) = spksGen;
end

 ttf_avg = mean(ttf_log,2);

 mov=squeeze(mean(mov,3));
 maskedMovdd= filterMov(mov,ones(size(mov,1),size(mov,2)),squeeze(ttf_avg));

  save('/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2008_08_27_6/On_large_type2.mat','cells','maskedMovdd','Y','ttf_log','ttf_avg','totalMaskAccept_log','-v7.3')
%% Test ASM one cell below ..



 interval=1;
%  idx = [1:end]
trainData=[1:floor(length(spksGen)*1)];
trainData_hr = [1:floor(length(spksGen_hr)*1)];

testData=[floor(length(spksGen)*0.8)+1 :length(spksGen)];
testData_hr=[floor(length(spksGen_hr)*0.8)+1 :length(spksGen_hr)];


 nSU=5;
 binnedResponsesbigd = spksGen(trainData);
binnedResponsesbigd_hr = spksGen_hr(trainData_hr);
 mov_use=maskedMovdd(:,trainData);
  filteredStimDim=size(mov_use,1); 
 
% [fitGMLM,output] = fitGMLM_MEL_EM(binnedResponsesbigd,mov_use,filteredStimDim,nSU,interval);  
 [fitGMLM,f_val(nSU)] = fitGMLM_MEL_EM_bias(binnedResponsesbigd,mov_use,filteredStimDim,nSU,interval); 
 fitGMLM_log{nSU} = fitGMLM;
 %[fitGMLM_2,output] = fitGMLM_MEL_EM_power2(binnedResponsesbigd,mov_use,filteredStimDim,nSU,interval,2);  
 % fitGMLM_log(ifit).fitGMLM = fitGMLM;  %%
%   [fitGMLM_full2,output]= fitGMLM_full(fitGMLM,binnedResponsesbigd_hr,mov_use);
%   fitGMLM_full2_log{nSU}=fitGMLM_full2;
  figure;
  plot(fitGMLM_full2.hist.hexpanded)

%% show filters

figure;

    
    fitGMLM_full2 = fitGMLM;
      mask = totalMaskAccept;
sta_dim1 = size(mask,1);
sta_dim2 = size(mask,2);
indexedframe = reshape(1:sta_dim1*sta_dim2,[sta_dim1,sta_dim2]);
masked_frame = indexedframe(logical(mask));

u_spatial_log = zeros(40,40,nSU);

for ifilt=1:nSU

u_spatial = reshape_vector(fitGMLM.Linear.filter{ifilt}(1:length(masked_frame)),masked_frame,indexedframe);
%u_spatial = reshape_vector(WN_uSq{ifilt},masked_frame,indexedframe);

subplot(3,2,ifilt)
imagesc(u_spatial(x_coord,y_coord));
colormap gray
colorbar
title(sprintf('GMLM Filter: %d',ifilt));
axis square
u_spatial = u_spatial.*totalMaskAccept;
u_spatial_log(:,:,ifilt) = u_spatial;
end
