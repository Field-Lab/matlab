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
WN_datafile = '2015-03-09-2/streamed/data038/data038';
WN_datafile_short='2015-03-09-2/streamed/data038/data038';
movie_xml = 'BW-8-2-0.48-11111-40x40';
stim_length=1800;% in seconds

%% Load stimuli


datarun=load_data(WN_datafile)
datarun=load_params(datarun)
datarun=load_sta(datarun);
datarun=load_neurons(datarun);

for cellID = 2371%datarun.cell_types{2}.cell_ids;
    cellID
%cell_glm_fit = sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_03_09_2/data038/CellType_OFF parasol/CellID_%d.mat',cellID);
%load(cell_glm_fit);
%%
extract_movie_response;

 %% EM like Max Expected Likelihood .. 
 interval=1;
%  idx = [1:end]
 binnedResponsesbigd = spksGen;
 mov_use=maskedMovdd;
 filteredStimDim=size(mov_use,1);
 nSU = 3;
 for ifit=1:50
 %[fitGMLM,output] = fitGMLM_MEL_EM(binnedResponsesbigd,mov_use,filteredStimDim,nSU,interval);  
 [fitGMLM,output] = fitGMLM_MEL_EM_bias(binnedResponsesbigd,mov_use,filteredStimDim,nSU,interval);  
 %[fitGMLM,output] = fitGMLM_MEL_EM_power2(binnedResponsesbigd,mov_use,filteredStimDim,nSU,interval,2);  
 fitGMLM_log(ifit).fitGMLM = fitGMLM;
 end
% save(sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2015_03_09_2/Cell%d',cellID),'fitGMLM_log','mov_use','binnedResponsesbigd','nSU','filteredStimDim','interval','totalMaskAccept2','totalMaskAccept','x_coord','y_coord');
 
 end
 
  %%
 [fitGMLM_full2,output]= fitGMLM_full(fitGMLM,spike.home,mov_use);
figure;
 plot(fitGMLM_full2.hist.hexpanded)
  %% Show learned filters;
  mask = totalMaskAccept;
sta_dim1 = size(mask,1);
sta_dim2 = size(mask,2);
indexedframe = reshape(1:sta_dim1*sta_dim2,[sta_dim1,sta_dim2]);
masked_frame = indexedframe(logical(mask));

u_spatial_log = zeros(40,40,nSU);

figure;
for ifilt=1:nSU
subplot(2,2,ifilt)
u_spatial = reshape_vector(fitGMLM.Linear.filter{ifilt}(1:length(masked_frame)),masked_frame,indexedframe);
imagesc(u_spatial(x_coord,y_coord));
colormap gray
colorbar
title(sprintf('GMLM Filter: %d',ifilt));
axis square
u_spatial = u_spatial.*totalMaskAccept;
u_spatial_log(:,:,ifilt) = u_spatial;
end

[v,I] = max(u_spatial_log,[],3);

xx=I.*(v>0.2);
xx=xx(x_coord,y_coord);
figure;
imagesc(xx);

%iso_response_bias_gmlm(binnedResponsesbigd,mov_use,fitGMLM);

%% Response prediction 

%% Load different movies
nConditions=6;
% make movies
interval=2;
condMov=cell(nConditions,1);
rawMovFrames=1272/(2);
icnt=0;
% make pixel histogram
for imov=[1,2,4,6,8,10]
    [stim,height,width,header_size] = get_raw_movie(sprintf('/Volumes/Data/2015-03-09-2/Visual/pc2015_03_09_2_data038/%d.rawMovie',imov),rawMovFrames,1);
    subtract_movies{3}=mean(stim,1);
    subtract_movies{3}=mean(stim,1)*0+127.5;
    movie=stim-repmat(subtract_movies{3},[rawMovFrames,1,1]);
    movie=movie/255;
    
    icnt=icnt+1;
    qq=permute(movie,[2,3,1]);
    ifcnt = 0;
    condMov{icnt}=zeros(size(qq,1),size(qq,2),size(qq,3)*interval);
    for iframe=1:size(qq,3)
        for irepeat=1:interval
            ifcnt=ifcnt+1;
            condMov{icnt}(:,:,ifcnt)=double(qq(:,:,iframe)); % cond mov is between 0 and 1 now!
        end
        
    end
    
end

%% Do predictions - full
% 
% % Load recorded response
% Null_datafile = '/Volumes/Analysis/2015-03-09-2/data042-from-data038_streamed_nps';
% neuronPath = [Null_datafile,sprintf('/data042-from-data038_streamed_nps.neurons')];
% condDuration=10.6;
% nConditions=6;
% cond_str=[];
% [spkColl,spkCondColl,h]=plot_raster_script_pc2015_03_09_2_light(cellID,nConditions,condDuration,cond_str,neuronPath);
% 
% 
% % Make predictions
% pred=cell(nConditions,1);
% for  icond=1:nConditions
% movd = condMov{icond};
%  maskedMov= filterMov(movd,mask,squeeze(tf));
%  maskedMov2=[maskedMov;ones(1,size(maskedMov,2))];
% nTrials=50;
%  pred{icond}= predictGMLM_full(fitGMLM_full2,maskedMov2,nTrials)';
% end
% 
% plot_record_prediction(spkCondColl,pred)

%% Do predictions - no feedback

% Load recorded response
Null_datafile = '/Volumes/Analysis/2015-03-09-2/data042-from-data038_streamed_nps';
neuronPath = [Null_datafile,sprintf('/data042-from-data038_streamed_nps.neurons')];
condDuration=10.6;
nConditions=6;
cond_str=[];
[spkColl,spkCondColl,h]=plot_raster_script_pc2015_03_09_2_light(cellID,nConditions,condDuration,cond_str,neuronPath);


% Make predictions
pred=cell(nConditions,1);
for  icond=1:nConditions
movd = condMov{icond};
 maskedMov= filterMov(movd,mask,squeeze(tf));
 maskedMov2=[maskedMov;ones(1,size(maskedMov,2))];
nTrials=50;
 pred{icond}= predictGMLM_bias(fitGMLM2,maskedMov,nTrials)';
  %pred{icond}= predictGMLM(fitGMLM,maskedMov,nTrials)';
 % pred{icond}= predictGMLM_gamma2(fitGMLM,maskedMov,nTrials,2)';
end

plot_record_prediction(spkCondColl,pred)