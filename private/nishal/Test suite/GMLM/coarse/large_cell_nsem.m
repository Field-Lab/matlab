
%% Verify cell ID in first and second sections match!
%% Find relevant pixels and tf using WN movie

% cellID=901;
% WN_datafile = '2012-09-27-3/data003/data003';
% movie_xml = 'RGB-10-2-0.48-11111';
% wrong_xml = 'RGB-10-2-0.48-22222';
% stim_length=1799;% 

cellID=1471;
WN_datafile = '2012-08-09-3/data002/data002';
movie_xml = 'RGB-8-1-0.48-11111';
wrong_xml = 'RGB-8-1-0.48-22222';
stim_length=1800;% 

datarun=load_data(WN_datafile)
datarun=load_params(datarun)
datarun=load_sta(datarun);
datarun=load_neurons(datarun);


 threshold_wrongSTA=sta_thr_using_wrong_movie(wrong_xml,datarun,stim_length,cellID);
 cell_params.thres = threshold_wrongSTA;

extract_movie_response3;

 mov_use_WN = maskedMovdd;
 binnedResponsesbigd_WN = spksGen;
filteredStimDim=size(mov_use,1);
nSU=9;
interval=1;
tic;
%[fitGMLM_WN,f_val] = fitGMLM_EM_bias(binnedResponsesbigd,mov_use,filteredStimDim,nSU,interval);
 [fitGMLM_WN,output] = fitGMLM_EM_power2(binnedResponsesbigd_WN,mov_use_WN,filteredStimDim,nSU,interval,2);  
toc;
   

      %% %% Load Large cell data for NSEM stimuli
      
%load('/Volumes/Lab/Users/akheitman/NSEM_Home/BlockedSpikes/2012-09-27-3/NSEM_mapPRJ/organizedspikes_Unknown_6902.mat');

% load('/Volumes/Lab/Users/akheitman/NSEM_Home/BlockedSpikes/2012-09-27-3/NSEM_mapPRJ/organizedspikes_OFFPar_901.mat');
% load('/Volumes/Lab/Users/akheitman/NSEM_Home/Stimuli/NSEM_eye-120-3_0-3600/fitmovie_schemeA_8pix_Identity_8pix.mat');

load('/Volumes/Lab/Users/akheitman/NSEM_Home/BlockedSpikes/2012-08-09-3/NSEM_mapPRJ/organizedspikes_OFFPar_1471.mat');
load('/Volumes/Lab/Users/akheitman/NSEM_Home/Stimuli/NSEM_eye-120-3_0-3600/fitmovie_schemeA_8pix_Identity_8pix.mat');

movie_frames = NSEMmovie.fitmovie.movie_byblock;

nmovies=length(movie_frames);
movie_full=zeros(80,40,7200*nmovies);
for imov=1:nmovies
    imov
movie_full(:,:,(imov-1)*7200+1:imov*7200) = double(movie_frames{imov}.matrix)/255 - 0.5;
end

spikes_by_block = organizedspikes.block.t_sp_withinblock(2:2:end);

binnedSpikeResponses =[];
for imov=1:nmovies
    
    ss{1}=spikes_by_block{imov}*20000; % Remove spikes which are too soon after change of movie.
    ss{1}(ss{1}<20000*30/120)=0;
    spkMat = makeSpikeMat(ss,1/120,7200);
    binnedSpikeResponses = [binnedSpikeResponses,spkMat];

end

%% Calculate STA

times = 1:length(binnedSpikeResponses);
t_sp=times(binnedSpikeResponses~=0);

STA=zeros(80,40,30);
for it_sp=t_sp(t_sp>40)
STA = STA+movie_full(:,:,it_sp:-1:it_sp-30+1);
end

STA =STA/numel(t_sp);

average_mov = repmat(mean(movie_full,3),[1,1,30]);


STA2 = STA - average_mov;
STA3 = STA./average_mov;

figure;
subplot(4,1,1);
imagesc(STA(:,:,6)');
colormap gray
axis image

subplot(4,1,2);
imagesc(STA2(:,:,6)');
colormap gray;
axis image

subplot(4,1,3);
imagesc(STA3(:,:,6)');
colormap gray;
axis image


subplot(4,1,4);
imagesc(-STA3(:,:,6)');
colormap gray;
axis image

%% Find relevant pixels and tf using NSEM
        stas{1}=zeros(80,40,3,30);
        stas{1}(:,:,1,:) = STA2;
        stas{1}(:,:,2,:) = STA2;
        stas{1}(:,:,3,:) = STA2;
        
         cell_params.STAlen=30;
         cell_params.thres=1.5;
       %[stas_clipped,totalMaskAccept2,CellMasks]= clipSTAs(stas,cell_params);
        [stas_clipped,totalMaskAccept2,CellMasks]= clipSTAs_largestblob(stas,cell_params);
        % fill in the bad entries in totalMaskAccept2
       
%         
%       ttf=squeeze(mean(mean(STA2.*repmat(totalMaskAccept2,[1,1,30]),1),2));
%       ttf=1000*ttf;
%       %tf=tf.*double(idx<15)';
%      
%       figure;
%       plot(ttf);
        
%% Filter with mask and tf
      mov = movie_full; 
      maskedMovdd= filterMov(mov,totalMaskAccept2,squeeze(ttf));
      
      %maskedMov2dd=[maskedMovdd;ones(1,size(maskedMovdd,2))];
      totalMaskAccept = totalMaskAccept2;


      %% Make fits
      mov_use = maskedMovdd;
      binnedResponsesbigd = binnedSpikeResponses';
      filteredStimDim = size(mov_use,1);
      interval=1;
      nSU=4;
      
      for ifit=1:1
          ifit
%      [fitGMLM,f_val(nSU)] = fitGMLM_MEL_EM_bias(binnedResponsesbigd,mov_use,filteredStimDim,nSU,interval); 
%      fitGMLM_log{ifit} = fitGMLM;
     
      [fitGMLM_NSEM,f_val] = fitGMLM_EM_bias(binnedResponsesbigd,mov_use,filteredStimDim,nSU,interval)
%  %[fitGMLM_2,output] = fitGMLM_MEL_EM_power2(binnedResponsesbigd,mov_use,filteredStimDim,nSU,interval,2);  
%  % fitGMLM_log(ifit).fitGMLM = fitGMLM;  %%
%   [fitGMLM_full2,output]= fitGMLM_full(fitGMLM,binnedResponsesbigd_hr,mov_use);
%   fitGMLM_full2_log{nSU}=fitGMLM_full2;
      end

      
      %% set GMLM
      
fitGMLM=fitGMLM_WN;

%% Analyse fit
mask = totalMaskAccept2;
sta_dim1 = size(mask,1);
sta_dim2 = size(mask,2);
indexedframe = reshape(1:sta_dim1*sta_dim2,[sta_dim1,sta_dim2]);
masked_frame = indexedframe(logical(mask));

u_spatial_log = zeros(sta_dim1,sta_dim2,nSU);
su_log = zeros(length(masked_frame));
figure;
for ifit=1:1
  %  fitGMLM=fitGMLM_log{ifit};
W=zeros(length(masked_frame),nSU);
for ifilt=1:nSU

u_spatial = reshape_vector(fitGMLM.Linear.filter{ifilt}(1:length(masked_frame)),masked_frame,indexedframe);

subplot(ceil((nSU+1)/2),ceil((nSU+1)/2),ifilt)
imagesc(u_spatial);
colormap gray
colorbar
title(sprintf('gmlm Filter: %d',ifilt));
axis image

u_spatial = u_spatial.*totalMaskAccept;
u_spatial_log(:,:,ifilt) = u_spatial;
W(:,ifilt) = fitGMLM.Linear.filter{ifilt}(1:length(masked_frame)); 
end


% make occurance histogram
 k_est=W';
 su=double((k_est)==repmat(max((k_est),[],1),[nSU,1]));
 su_log = su_log + su'*su;
end


%% Response prediction

spikes_by_block_test = organizedspikes.block.t_sp_withinblock(1:2:end);

testSpikeResponses =[];
for imov=1:length(spikes_by_block_test)
    
    ss{1}=spikes_by_block_test{imov}*20000; % Remove spikes which are too soon after change of movie.
    ss{1}(ss{1}<20000*30/120)=0;
    spkMat = makeSpikeMat(ss,1/120,3600);
    testSpikeResponses = [testSpikeResponses;spkMat];
end

figure;
[x1,y1]=plotSpikeRaster(testSpikeResponses~=0,'PlotType','vertline');

% test movie? 
test=load('/Volumes/Lab/Users/akheitman/NSEM_Home/Stimuli/NSEM_eye-120-3_0-3600/testmovie_schemeA_8pix_Identity_8pix.mat');
test_mov = double(test.testmovie.matrix)/255 - 0.5;

% Do same processing as it was done to stimuli while fitting it.
mov = test_mov; 

if(size(mov,1)~=size(totalMaskAccept2,1))
mov2 = zeros(size(totalMaskAccept2,1),size(totalMaskAccept2,2),size(mov,3));

for iframe=1:size(mov,3)
mov2(:,:,iframe) = imresize(mov(:,:,iframe),size(totalMaskAccept2,1)/size(mov,1));
end
mov=mov2;
end

maskedMov_test= filterMov(mov,totalMaskAccept2,squeeze(ttf));


% make response prediction.
nTrials=59;
interval=1;
%[pred,lam]= predictGMLM_bias(fitGMLM,maskedMov_test,nTrials,interval);
[pred,lam] = predictGMLM_gamma2(fitGMLM,maskedMov_test,nTrials,2,interval);

pred=pred';
figure;
[x2,y2] = plotSpikeRaster(pred~=0,'plotType','vertline');

figure;
plot(x1*30/max(x1(:)),y1+max(y2(:)),'r');
hold on;
plot(x2*30/max(x2(:)),y2,'k');
