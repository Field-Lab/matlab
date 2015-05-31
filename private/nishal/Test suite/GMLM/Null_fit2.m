%% Null movie fit GMLM

% Dataset 2015-02-24-5

%% data000 and get tf and important pixels

WN_datafile = '2015-02-24-5/streamed/data006/data006';
movie_xml = 'BW-20-8-0.48-11111-16x16';
stim_length=900;% in seconds
cellID =286;

datarun=load_data(WN_datafile)
datarun=load_params(datarun)
datarun=load_sta(datarun);
datarun=load_neurons(datarun);

extract_movie_response2;

%% analyze data004,data000
% Condition strings
nConditions=2;
cond_str=cell(3,1);
cond_str{1}='Original';
cond_str{2}='Cell group 1 Spatial';
cond_str{3}='OFF parasol';
interestingConditions=[1,2,3];



%% Load spikes
WN_datafile = '2015-02-24-5/streamed/data006/data006';
Null_datafile = '/Volumes/Analysis/2015-02-24-5/data011-from-data006_s_nps';
%WN_mov='/Volumes/Analysis/stimuli/white-noise-xml/RGB-10-2-0.48-11111.xml';
neuronPath = [Null_datafile,sprintf('/data011-from-data006_s_nps.neurons')];


datarun=load_data(WN_datafile)
datarun=load_params(datarun)

condDuration=92;
nConditions=2;

    % [spkColl,spkCondColl,h]=plot_raster_script_pc2015_02_24_2(datarun,WN_datafile,WN_datafile_full,Null_datafile,InterestingCell_vis_id,imov,ref_cell_number,nConditions,condDuration,cond_str,neuronPath);
% [spkColl,spkCondColl,h]=plot_raster_script_pc2015_03_09_3_light(cellID,nConditions,condDuration,cond_str,neuronPath);


 [spkColl,spkCondColl,h]=plot_raster_script_pc2015_02_24_2_light(cellID,nConditions,condDuration,cond_str,neuronPath)
%% data004 movies
% make movies
interval=8;
condMov=cell(nConditions,1);
rawMovFrames=22080/(interval);
icnt=0;
% make pixel histogram
for imov=[20]
    [stim,height,width,header_size] = get_raw_movie(sprintf('/Volumes/Data/2015-02-24-5/Visual/pc2015_02_24_5_data006/%d.rawMovie',imov),rawMovFrames,1);
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
            condMov{icnt}(:,:,ifcnt)=qq(:,:,iframe); % cond mov is between 0 and 1 now!
        end
        
    end
    
end
movq=condMov{1};
condMov{1}=movq(:,:,1:floor(end/2));
condMov{2}=movq(:,:,floor(end/2)+1:end);

%%
useCond=2;
spksGen =  makeSpikeMat(spkCondColl(useCond).spksColl,1/120,condDuration*120);
ss=spksGen';
spksGen1 = ss(:);
movNull1=repmat(condMov{useCond},[1,1,5]);
spksGen1(size(movNull1,3))=0;

useCond=1;
spksGen =  makeSpikeMat(spkCondColl(useCond).spksColl,1/120,condDuration*120);
ss=spksGen';
spksGen2 = ss(:);
movNull2=repmat(condMov{useCond},[1,1,5]);
spksGen2(size(movNull2,3))=0;

spksGen=[spksGen1;spksGen2];
movNull = zeros(size(movNull1,1),size(movNull1,2),2*size(movNull1,3));
movNull(:,:,1:end/2) = movNull1;
movNull(:,:,(end/2)+1:end)=movNull2;
%% Fit 

 maskedMovdd= filterMov(movNull,totalMaskAccept2,squeeze(tf));
 
  interval=1;
%  idx = [1:end]
trainData=[1:floor(length(spksGen)*1)];
%testData=[floor(length(spksGen)*0.8)+1 :length(spksGen)];
 binnedResponsesbigd = spksGen(trainData);
 mov_use=maskedMovdd(:,trainData);
 
 filteredStimDim=size(mov_use,1); 

 for nSU =1;
 binnedResponsesbigd = spksGen(trainData);
 mov_use=maskedMovdd(:,trainData);
 
 for ifit=1:1
% [fitGMLM,output] = fitGMLM_MEL_EM(binnedResponsesbigd,mov_use,filteredStimDim,nSU,interval);  
 [fitGMLM,f_val(nSU)] = fitGMLM_MEL_EM_bias(binnedResponsesbigd,mov_use,filteredStimDim,nSU,interval);  
 fitGMLM_sulog{nSU}=fitGMLM;
 %[fitGMLM,output] = fitGMLM_MEL_EM_power2(binnedResponsesbigd,mov_use,filteredStimDim,nSU,interval,2);  
 % fitGMLM_log(ifit).fitGMLM = fitGMLM;  %%
 [fitGMLM_full2,output]= fitGMLM_full(fitGMLM,binnedResponsesbigd,mov_use);
 %fitGMLM_full2_log{nSU}=fitGMLM_full2;
figure;
plot(fitGMLM_full2.hist.hexpanded)
 end
 
 end
% save('/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2015_03_09_3/fits1939','fitGMLM_sulog','fitGMLM_full2_log');
%  % Training Err 
%   nTrials=1;
%   binnedResponsesbigd_train = spksGen(trainData);
%  mov_use_train=maskedMovdd(:,trainData);
%  [pred_train{nSU},lam_train{nSU}]= predictGMLM_bias(fitGMLM,mov_use_train,nTrials,interval); 
%  rec=binnedResponsesbigd_train;
%  [f_val_train(nSU),R2_train(nSU)]=f_r2_test(lam_train{nSU},rec,interval);
%  
%  % Testing.
%  nTrials=1;
%   binnedResponsesbigd_test = spksGen(testData);
%  mov_use_test=maskedMovdd(:,testData);
%  [pred_test{nSU},lam_test{nSU}]= predictGMLM_bias(fitGMLM,mov_use_test,nTrials,interval); 
%  rec=binnedResponsesbigd_test;
%  [f_val_test(nSU),R2_test(nSU)]=f_r2_test(lam_test{nSU},rec,interval);
% 
%  end
% % save(sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2015_03_09_2/Cell%d_full',cellID),'fitGMLM_log','fitGMLM_full2','mov_use','binnedResponsesbigd','nSU','filteredStimDim','interval','totalMaskAccept2','totalMaskAccept','x_coord','y_coord');
% 
% 
% figure('Color','w');
% plot(1:length(f_val_train),f_val_train,'--*');
% hold on;
% plot(1:length(f_val_test),f_val_test,'--*');
% legend('Training','Testing');
% xlabel('number of SU');
% ylabel('negative Likelihood');
% title('Training and Testing Likelihood value');
% 
% figure('Color','w');
% plot(1:length(R2_train),R2_train,'--*');
% hold on;
% plot(1:length(R2_test),R2_test,'--*');
% legend('Training','Testing');
% xlabel('number of SU');
% ylabel('R2');
% title('Training and Testing R-sq value');

%% 
nSU=2;
%fitGMLM_full2=fitGMLM_full2_log{useSU};
fitGMLM=fitGMLM_sulog{nSU};

 mask = totalMaskAccept;
sta_dim1 = size(mask,1);
sta_dim2 = size(mask,2);
indexedframe = reshape(1:sta_dim1*sta_dim2,[sta_dim1,sta_dim2]);
masked_frame = indexedframe(logical(mask));

u_spatial_log = zeros(dim1,dim2,nSU);

figure;
for ifilt=1:nSU
subplot(2,2,ifilt)
u_spatial = reshape_vector(fitGMLM_full2.Linear.filter{ifilt}(1:length(masked_frame)),masked_frame,indexedframe);
imagesc(u_spatial(x_coord,y_coord));
colormap gray
colorbar
title(sprintf('GMLM Filter: %d',ifilt));
axis square
u_spatial = u_spatial.*totalMaskAccept;
u_spatial_log(:,:,ifilt) = u_spatial;
end

[v,I] = max(u_spatial_log,[],3);

xx=I.*(v>0.1);
xx=xx(x_coord,y_coord);
figure;
imagesc(xx);

%iso_response_bias_gmlm(binnedResponsesbigd,mov_use,fitGMLM);

%su_activation_plot(fitGMLM_full2,mov_use);
%% data003 movies
% make movies
interval=8;
condMov2=cell(nConditions,1);
rawMovFrames=4320/(interval);
icnt=0;
% make pixel histogram
for imov=[19]
    [stim,height,width,header_size] = get_raw_movie(sprintf('/Volumes/Data/2015-02-24-5/Visual/pc2015_02_24_5_data006/%d.rawMovie',imov),rawMovFrames,1);
    subtract_movies{3}=mean(stim,1);
    subtract_movies{3}=mean(stim,1)*0+127.5;
    movie=stim-repmat(subtract_movies{3},[rawMovFrames,1,1]);
    movie=movie/255;
    
    icnt=icnt+1;
    qq=permute(movie,[2,3,1]);
    ifcnt = 0;
    condMov2{icnt}=zeros(size(qq,1),size(qq,2),size(qq,3)*interval);
    for iframe=1:size(qq,3)
        for irepeat=1:interval
            ifcnt=ifcnt+1;
            condMov2{icnt}(:,:,ifcnt)=qq(:,:,iframe); % cond mov is between 0 and 1 now!
        end
        
    end
    
end
movqq=condMov2{1};

condMov2{1}=movqq(:,:,1:floor(end/3));
condMov2{2}=movqq(:,:,floor(end/3)+1:floor(2*end/3));
condMov2{3}=movqq(:,:,floor(2*end/3)+1:end);
%% data003 spikes
Null_datafile = '/Volumes/Analysis/2015-02-24-5/data010-from-data006_s_nps';
neuronPath = [Null_datafile,sprintf('/data010-from-data006_s_nps.neurons')];
nConditions=3;
condDuration=12;
 [spkColl,spkCondColl,h]=plot_raster_script_pc2015_02_24_2_light(cellID,nConditions,condDuration,cond_str,neuronPath)
%% make prediction

% Make predictions
pred=cell(nConditions,1); lam=cell(nConditions,1);
for  icond=1:nConditions
 maskedMov= filterMov(condMov2{icond},mask,squeeze(tf));
 %maskedMov2=[maskedMov;ones(1,size(maskedMov,2))];
nTrials=30;
interval=1;
 %[pred{icond},lam{icond}]= predictGMLM_bias(fitGMLM,maskedMov,nTrials,interval);
pred{icond} = predictGMLM_full(fitGMLM_full2,maskedMov,nTrials);
  %pred{icond}= predictGMLM(fitGMLM,maskedMov,nTrials)';
 % pred{icond}= predictGMLM_gamma2(fitGMLM,maskedMov,nTrials,2)';
  pred{icond}=pred{icond}';
end

plot_record_prediction(spkCondColl,pred)      

%plot_record_prediction3(spkCondColl,pred1cell,pred10cell);