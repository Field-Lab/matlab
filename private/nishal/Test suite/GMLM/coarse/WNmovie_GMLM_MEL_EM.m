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

for cellID = [datarun.cell_types{2}.cell_ids];
    cellID
%cell_glm_fit = sprintf('/Volumes/Lab/Users/bhaishahster/analyse_2015_03_09_2/data038/CellType_OFF parasol/CellID_%d.mat',cellID);
%load(cell_glm_fit);
%%
extract_movie_response2;

 %% EM like Max Expected Likelihood .. 
 interval=1;
%  idx = [1:end]
trainData=[1:floor(length(spksGen)*1)];
trainData_hr = [1:floor(length(spksGen_hr)*1)];

testData=[floor(length(spksGen)*0.8)+1 :length(spksGen)];
testData_hr=[floor(length(spksGen_hr)*0.8)+1 :length(spksGen_hr)];


 for nSU =4:5%1:2 %1:10%1:filteredStimDim
 binnedResponsesbigd = spksGen(trainData);
binnedResponsesbigd_hr = spksGen_hr(trainData_hr);
 mov_use=maskedMovdd(:,trainData);
  filteredStimDim=size(mov_use,1); 
 
  for ifit=1:1
% [fitGMLM,output] = fitGMLM_MEL_EM(binnedResponsesbigd,mov_use,filteredStimDim,nSU,interval);  
 [fitGMLM,f_val(nSU)] = fitGMLM_MEL_EM_bias(binnedResponsesbigd,mov_use,filteredStimDim,nSU,interval); 
 fitGMLM_log{nSU} = fitGMLM;
 %[fitGMLM_2,output] = fitGMLM_MEL_EM_power2(binnedResponsesbigd,mov_use,filteredStimDim,nSU,interval,2);  
 % fitGMLM_log(ifit).fitGMLM = fitGMLM;  %%
  [fitGMLM_full2,output]= fitGMLM_full(fitGMLM,binnedResponsesbigd_hr,mov_use);
  fitGMLM_full2_log{nSU}=fitGMLM_full2;
  figure;
  plot(fitGMLM_full2.hist.hexpanded)
 end
 end
 [WNSTA,WNSTC,WN_uSq]=compute_STA_STC(binnedResponsesbigd,mov_use);
 save(sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2015_03_09_2/data038/Off parasol/Cell%d_full_su4_5.mat',cellID),'fitGMLM_log','fitGMLM_full2_log','mov_use','binnedResponsesbigd','filteredStimDim','interval','totalMaskAccept2','totalMaskAccept','x_coord','y_coord','WNSTA','WNSTC','WN_uSq');
%% Compute STC 

end
%% Compare training and test set likelihoods accross models when true spikes known.

 for nSU=1:4
     nSU
     fitGMLM_full2=fitGMLM_full2_log{nSU};
 % Training Err 
  nTrials=100;
  binnedResponsesbigd_train = spksGen(trainData);
binnedResponsesbigd_hr_train = spksGen_hr(trainData_hr);
 mov_use_train=maskedMovdd(:,trainData);
  [nlikelihood,lam]= likelihood_fitGMLM_full(fitGMLM_full2,binnedResponsesbigd_hr_train,mov_use_train);
f_val_train(nSU)=nlikelihood;
lam_train{nSU}=lam;
 
 % Testing.
 nTrials=100;
  binnedResponsesbigd_test = spksGen(testData);
  binnedResponsesbigd_hr_test = spksGen_hr(testData_hr);
 mov_use_test=maskedMovdd(:,testData);
   [nlikelihood,lam]= likelihood_fitGMLM_full(fitGMLM_full2,binnedResponsesbigd_hr_test,mov_use_test);
 
f_val_test(nSU)=nlikelihood;
lam_test{nSU} = lam;
 end
 
 %% Compare training and test set likelihoods across models when true spikes NOT known.
 for nSU=1:3
     nSU
     fitGMLM_full2=fitGMLM_full2_log{nSU};
 % Training Err 
  nTrials=100;
  binnedResponsesbigd_train = spksGen(trainData);
binnedResponsesbigd_hr_train = spksGen_hr(trainData_hr);
 mov_use_train=maskedMovdd(:,trainData);
 %[pred,lam]= predictGMLM_bias(fitGMLM,mov_use_train,nTrials,interval); 
 [pred,lam] = predictGMLM_full(fitGMLM_full2,mov_use_train,nTrials);
% lam=lam*1200;
 lam = mean(pred,2)'*1200; lam(lam==0)=0.001;
 rec=binnedResponsesbigd_hr_train;
 [f_val_train(nSU),R2_train(nSU)]=f_r2_test(lam,rec,interval);
 lam_train{nSU}=lam;
 
 % Testing.
 nTrials=100;
  binnedResponsesbigd_test = spksGen(testData);
  binnedResponsesbigd_hr_test = spksGen_hr(testData_hr);
 mov_use_test=maskedMovdd(:,testData);
% [pred,lam]= predictGMLM_bias(fitGMLM,mov_use_test,nTrials,interval); 
  %[pred,lam] = predictGMLM_full(fitGMLM_full2,mov_use_test,nTrials);
  % [nlikelihood]= likelihood_fitGMLM_full(fitGMLM_full2,binnedResponsesbigd_hr_test,mov_use_test);
 % lam=lam*1200;
   lam = mean(pred,2)'*1200; lam(lam==0)=0.001;
 rec=binnedResponsesbigd_hr_test;
 [f_val_test(nSU),R2_test(nSU)]=f_r2_test(lam,rec,interval);
 lam_test{nSU}=lam;
 end
 %% Plot training and test set likelihood and R-2 values.
 
figure('Color','w');
plot(1:length(f_val_train),f_val_train,'--*');
hold on;
plot(1:length(f_val_test),f_val_test,'--*');
legend('Training','Testing');
xlabel('number of SU');
ylabel('negative Likelihood');
title('Training and Testing Likelihood value');

figure('Color','w');
plot(1:length(R2_train),R2_train,'--*');
hold on;
plot(1:length(R2_test),R2_test,'--*');
legend('Training','Testing');
xlabel('number of SU');
ylabel('R2');
title('Training and Testing R-sq value');

 
%% Prediction 
 
  %% Show learned filters;
figure;
for nSU=1:4
    
    fitGMLM_full2 = fitGMLM_full2_log{nSU};
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

[v,I] = max(u_spatial_log,[],3);

xx=I.*(v>0.2);
xx=xx(x_coord,y_coord);

subplot(2,2,nSU);
imagesc(repelem(xx,10,10));
title(sprintf('Num SU : %d',nSU));
hold on;
end

%iso_response_bias_gmlm(binnedResponsesbigd,mov_use,fitGMLM);

%su_activation_plot(fitGMLM,mov_use);
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

% Load recorded response
Null_datafile = '/Volumes/Analysis/2015-03-09-2/data042-from-data038_streamed_nps';
neuronPath = [Null_datafile,sprintf('/data042-from-data038_streamed_nps.neurons')];
condDuration=10.6;
nConditions=6;
cond_str=[];
[spkColl,spkCondColl,h]=plot_raster_script_pc2015_03_09_2_light(cellID,nConditions,condDuration,cond_str,neuronPath);


% Make predictions
pred1=cell(nConditions,1);
for  icond=1:nConditions
movd = condMov{icond};
 maskedMov= filterMov(movd,mask,squeeze(tf));
 maskedMov2=[maskedMov;ones(1,size(maskedMov,2))];
nTrials=30;
 pred1{icond}= predictGMLM_full(fitGMLM_full2_log{1},maskedMov,nTrials)';
end

% Make predictions
pred3=cell(nConditions,1);
for  icond=1:nConditions
movd = condMov{icond};
 maskedMov= filterMov(movd,mask,squeeze(tf));
 maskedMov2=[maskedMov;ones(1,size(maskedMov,2))];
nTrials=30;
 pred3{icond}= predictGMLM_full(fitGMLM_full2_log{3},maskedMov,nTrials)';
end

figure;
plot_record_prediction(spkCondColl,pred3);

figure;
plot_record_prediction3(spkCondColl,pred3,pred1);
%% Do predictions - no feedback

% Load recorded response
Null_datafile = '/Volumes/Analysis/2015-03-09-2/data042-from-data038_streamed_nps';
neuronPath = [Null_datafile,sprintf('/data042-from-data038_streamed_nps.neurons')];
condDuration=10.6;
nConditions=6;
cond_str=[];
[spkColl,spkCondColl,h]=plot_raster_script_pc2015_03_09_2_light(cellID,nConditions,condDuration,cond_str,neuronPath);


% Make predictions
pred=cell(nConditions,1); lam=cell(nConditions,1);
for  icond=1:nConditions
movd = condMov{icond};
 maskedMov= filterMov(movd,mask,squeeze(tf));
 maskedMov2=[maskedMov;ones(1,size(maskedMov,2))];
nTrials=50;
interval=1;
 [pred{icond},lam{icond}]= predictGMLM_bias(fitGMLM,maskedMov,nTrials,interval);
 pred{icond}=pred{icond}';
  %pred{icond}= predictGMLM(fitGMLM,maskedMov,nTrials)';
 % pred{icond}= predictGMLM_gamma2(fitGMLM,maskedMov,nTrials,2)';
end

plot_record_prediction(spkCondColl,pred)      


for icond=1:6
       predd=lam{icond}; 
     
       pred_ss = zeros(length(predd)/10,1);
       for itime=1:length(pred_ss)
       pred_ss(itime) = sum(predd((itime-1)*10+1:(itime)*10));
       end
       
       pred_ss=repmat(pred_ss,[29,1]);
       
       rec = makeSpikeMat(spkCondColl(icond).spksColl,1/120,condDuration/(1/120));
       rec=rec';
       rec=rec(:);
       
       
       % R2 value method 2
       x1 = pred_ss; y1 = rec; n=length(y1);
       r = (n*x1'*y1 - sum(x1)*sum(y1))/(sqrt(n*sum(x1.^2) - sum(x1)^2) * sqrt(n*sum(y1.^2) - sum(y1)^2));
       R2_log(icond) = r^2;
       
end
R2_log


%% compare model predictions on test data in two different conditions.

lam1=lam_test{1};
lam2=mean(lam1)*ones(1,length(lam1));%lam_test{3};
binnedResponsesbigd_hr_test = spksGen_hr(testData_hr);

[like_full,like_diff,diffll0]=compare_lam_likelihood(lam1,lam2,binnedResponsesbigd_hr_test);
figure;
[X0,N0]=hist(diffll0,1000);



lam1=lam_test{1};
lam2=lam_test{3};
binnedResponsesbigd_hr_test = spksGen_hr(testData_hr);

[like_full,like_diff,diffll]=compare_lam_likelihood(lam1,lam2,binnedResponsesbigd_hr_test);
hold on
[X,N]= hist(diffll,1000);

figure;
semilogy(N0,X0/sum(X0));
hold on;
semilogy(N,X/sum(X));
%legend(sprintf('GLM 1-Null %f',(diffll0<0)/length(diffll0)),sprintf('GLM3-GLM1 %f',(diffll<0)/length(diffll)))
