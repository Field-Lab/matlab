% Null space simulation suite ! 

%% Generate RFs - subunit model
nSubunits = 4;
Filtdim1 = 6;
Filtdim2 = 6;
Filtlen = 30;

subunits=cell(nSubunits,1);

k=zeros(Filtdim1,Filtdim2,1,Filtlen);
k(320-318,160-158,:,:)=1;
k(323-318,162-158,:,:)=1;
k(323-318,163-158,:,:)=1;
subunit_scale=1%1;
subunits{1}=k*subunit_scale;

k=zeros(Filtdim1,Filtdim2,1,Filtlen);
k(320-318,159-158,:,:)=1;
subunit_scale=1%1.2;
subunits{2}=k*subunit_scale;

k=zeros(Filtdim1,Filtdim2,1,Filtlen);
k(322-318,161-158,:,:)=1;
subunit_scale=1%0.5;
subunits{3}=k*subunit_scale;

k=zeros(Filtdim1,Filtdim2,1,Filtlen);
k(322-318,160-158,:,:)=1;
k(319-318,159-158,:,:)=1;
subunit_scale=1%1.5;
subunits{4}=k*subunit_scale;


figure;
for isubunit=1:nSubunits
subplot(2,2,isubunit);
imagesc(subunits{isubunit}(:,:,1,4));
colormap gray
title(sprintf('Subunit: %d',isubunit));
end

figure;
mask=double((subunits{1}(:,:,1,4)+subunits{2}(:,:,1,4)+subunits{3}(:,:,1,4)+subunits{4}(:,:,1,4))~=0);
mask(323-318,162-158,:,:)=1;
imagesc(mask);
% Temporal properties?

scale_one=1;
scale_two=0.25;
tau_one=4;
tau_two=10;
n_filters=6;
t=[0:29];
tf = scale_one*((t/tau_one).^n_filters).*exp(-n_filters*(t/tau_one -1)) - scale_two*((t/tau_two).^n_filters).*exp(-n_filters*(t/tau_two -1));
figure;
plot(tf)

tf2=zeros(1,1,1,Filtlen);
tf2(1,1,1,:)=tf;
tf=tf2;
clear tf2
tfRep=repmat(tf,[Filtdim1,Filtdim2,1,1]);
title('Temporal Filter');

subunits{1}=subunits{1}.*tfRep;
subunits{2}=subunits{2}.*tfRep;
subunits{3}=subunits{3}.*tfRep;
subunits{4}=subunits{4}.*tfRep;


% sub-unit weights
subunitWeights=zeros(nSubunits,1);
subunitWeights(1)=1%1;
subunitWeights(2)=1%1.2;
subunitWeights(3)=1%0.7;
subunitWeights(4)=1%1.5;

% sub-unit non-linearity

f= @(x) double(x>0).*(1.6*x);
%f=@(x) exp(0.6*x);
% Ganglion cell non-linearity
%N= @(x) exp(0.15*(x));
N=@(x) double(x>0)*(0.02).*(3.4*x).^2 ; % use it! , 2 is DC current
%N=@(x) x-min(x(:));
%N = @(x) 15./(1+exp(-1.5*(x-5)));
% 
% 
% figure;
% for itime=1:30
%     itime
% for isubunit=1:nSubunits
% subplot(2,2,isubunit);
% imagesc(subunits{isubunit}(:,:,1,itime));
% caxis([min(subunits{isubunit}(:)),max(subunits{isubunit}(:))]);
% colormap gray
% colorbar
% %title(sprintf('Subunit: %d',isubunit));
% 
% end
% % pause
% end

%% Generate white noise movie
 movieLen=120*30*60;
mov=zeros(Filtdim1,Filtdim2,movieLen);
movie_idx=2; 
if(movie_idx==1)
mov(319-310,158-150,:)=0.5; % Stimulate a sub-unit
mov(320-310,160-150,:)=0.5;
mov=mov+rand(Filtdim1,Filtdim2,movieLen)*0.3;
end

if(movie_idx==2)
mov=double(rand(Filtdim1,Filtdim2,movieLen)>0.5)-0.5;
end

mov(:,:,1:30)=0;

figure;
for itime=40:50
   
imagesc(mov(:,:,itime));
colormap gray
colorbar
axis image
caxis([-0.5,0.5]);
pause(0.01)
end

%% Generate responses
% Calculate filter output for each sub-unit for each frame and calculate
% number of spikes for each frame-bin (binned response) .. So that would be
% used for STA calculation ? 



mov2=zeros(Filtdim1 ,Filtdim2,movieLen+Filtlen-1);
mov2(:,:,Filtlen:movieLen+Filtlen-1)=mov; % Append zeros before the movie
nTrials=1;
SubUnit_Response_test_movie_script

% Spike triggered sub-unit Input 
spikeTriggeredSubUnitInput(binnedResponses,cell_resp)

 %% Calculate STA


% My own STA code 
STA=zeros(Filtdim1,Filtdim2,Filtlen);
useTrial=1;
for iframe=30:movieLen
STA=STA+mov(:,:,iframe:-1:iframe-Filtlen+1)*binnedResponses(iframe(:,useTrial));
end
STA=STA/sum(binnedResponses(:,useTrial));


figure
 for itime=[1:Filtlen]
 imagesc(squeeze(STA(:,:,itime)));colormap gray
 caxis([min(STA(:)),max(STA(:))]);
 colorbar
 pause(1/120)
 end
 
%% Fit GMLM_afterSTC1 - detailed version .. 
%fitGMLM = fitGMLM_afterSTC(binnedResponses,mov,WN_uSq,WNSTA,subunits);

%% Train GMLM..
 maskedMov= filterMov(mov,mask,squeeze(tf));
 maskedMov2=[maskedMov;ones(1,size(maskedMov,2))];
% [fitGMLM,output] = fitGMLM_afterSTC_simplified(binnedResponses,maskedMov,7,4);
 
 %% EM like Max Expected Likelihood .. 
 interval=1;
 [fitGMLM,output] = fitGMLM_MEL_EM(binnedResponses,maskedMov2,8,4,interval);   
 idx=1:length(binnedResponses);
 spk_time = idx(binnedResponses>0)/120;
 [fitGMLM_full,output]= fitGMLM_full(fitGMLM,spk_time,maskedMov);
 %% Show learned filters;
sta_dim1 = size(mask,1);
sta_dim2 = size(mask,2);
indexedframe = reshape(1:sta_dim1*sta_dim2,[sta_dim1,sta_dim2]);
masked_frame = indexedframe(logical(mask));

figure;
for ifilt=1:4
subplot(2,2,ifilt)
u_spatial = -reshape_vector(fitGMLM.Linear.filter{ifilt}(1:length(masked_frame)),masked_frame,indexedframe);
imagesc(u_spatial);
colormap gray
colorbar
title(sprintf('GMLM Filter: %d',ifilt));
end

 
%% Generate white noise movie
 movieLen=120*10;
mov=zeros(Filtdim1,Filtdim2,movieLen);
movie_idx=2; 
if(movie_idx==1)
mov(319-310,158-150,:)=0.5; % Stimulate a sub-unit
mov(320-310,160-150,:)=0.5;
mov=mov+rand(Filtdim1,Filtdim2,movieLen)*0.3;
end

if(movie_idx==2)
mov=double(rand(Filtdim1,Filtdim2,movieLen)>0.5)-0.5;
end

mov(:,:,1:30)=0;

figure;
for itime=40:50
   
imagesc(mov(:,:,itime));
colormap gray
colorbar
axis image
caxis([-0.5,0.5]);
pause(0.01)
end

%% Generate responses
% Calculate filter output for each sub-unit for each frame and calculate
% number of spikes for each frame-bin (binned response) .. So that would be
% used for STA calculation ? 



mov2=zeros(Filtdim1 ,Filtdim2,movieLen+Filtlen-1);
mov2(:,:,Filtlen:movieLen+Filtlen-1)=mov; % Append zeros before the movie
nTrials=30;
SubUnit_Response_test_movie_script

%% Predict response using model
 maskedMov= filterMov(mov,mask,squeeze(tf));
 maskedMov2=[maskedMov;ones(1,size(maskedMov,2))];
nTrials=30;
 predictedResponses = predictGMLM(fitGMLM,maskedMov,nTrials);
 
 % Compare actual and predicted
figure;
[x1,y1]=plotSpikeRaster(binnedResponses'>0,'PlotType','vertline');

[x2,y2]=plotSpikeRaster(predictedResponses'>0,'PlotType','vertline');
figure;
plot(x1,y1+max(y2),'k');
hold on;
plot(x2/10,y2,'r');
legend('Recorded','Predicted');

%% Nulling experiment with GMLM
movieLen=120*15;
null_compute_subUnit_test
 
nTrials=50;
analyse_null_subUnit_ts

% Compare actual and predicted for original movie.
mov=movOrig(:,:,Filtlen:movie_new_len+Filtlen-1);
binnedResponses = binnedResponseOrig;
 maskedMov= filterMov(mov,mask,squeeze(tf));
 maskedMov2=[maskedMov;ones(1,size(maskedMov,2))];
nTrials=50;
 predictedResponses = predictGMLM(fitGMLM,maskedMov2,nTrials);
 

[x1,y1]=plotSpikeRaster(binnedResponses'>0,'PlotType','vertline');

[x2,y2]=plotSpikeRaster(predictedResponses'>0,'PlotType','vertline');

 % Compare actual and predicted for Null movie.
mov=movNull(:,:,Filtlen:movie_new_len+Filtlen-1);
binnedResponses = binnedResponseNull;
 maskedMov= filterMov(mov,mask,squeeze(tf));
 maskedMov2=[maskedMov;ones(1,size(maskedMov,2))];
nTrials=50;
 predictedResponses = predictGMLM(fitGMLM,maskedMov2,nTrials);
 
[x3,y3]=plotSpikeRaster(binnedResponses'>0,'PlotType','vertline');

[x4,y4]=plotSpikeRaster(predictedResponses'>0,'PlotType','vertline');

 % Compare actual and predicted for Original movie

figure;
plot(x1,y1+3*max(y2),'r');
hold on;
plot(x2/10,y2+2*max(y2),'k');
hold on;

hold on;
plot(x3,y3+max(y2),'r');
hold on;
plot(x4/10,y4,'k');
ylim([0,4*max(y2)])
xlim([0,max(x1)]);
title('Original & Null movie');
set(gca,'ytick',[]);
set(gca,'xtick',[]);
legend('Recorded Original ','Predicted Original','Recorded Null','Predicted Null');

