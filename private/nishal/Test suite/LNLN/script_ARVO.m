
%% Make model 
model = model_LNLN();

%% Generate test response

pixelSz=16;
sz = model.gridSzX/ (pixelSz*3);

% 
% Tlen = 120*10;
% movie = (randn(sz,sz,Tlen)>0)-0.5;
% dt=1/120;
% nTrials=30;
% [response,~] = generateResp_LNLN(model,movie,dt,nTrials);
% h= figure;
% plotSpikeRaster(response~=0,'PlotType','vertline');
% %print(h,'-dpdf',sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/modelCellLNLN/model1/sample_firing.pdf'));

Tlen = 120*30*30;
movie = (randn(sz,sz,Tlen)>0)-0.5;
dt=1/120;
nTrials=1;
[response,~] = generateResp_LNLN(model,movie,dt,nTrials);
%% Calulate STA 
idx = 1:Tlen;
spktm = idx(response~=0 & idx >30);

STA= zeros(sz,sz,30);


for itime=spktm
STA = STA + movie(:,:,itime-29:itime);
end

STA = STA / length(spktm);

figure;
for itime=26
    itime
imagesc(STA(:,:,itime));
colormap gray
caxis([min(STA(:)),max(STA(:))]);
%pause
end

% STA strongest frame with cone map;
strongestFrame = STA(:,:,26);
szstr = size(strongestFrame,1);

ssf = repelem(strongestFrame,model.gridSzX/szstr,model.gridSzX/szstr,3);

h=figure;
imagesc( repelem(sum(model.totalConeMap3D,3)==0,1,1,3).*ssf*30 + model.totalConeMap3D);
axis image
title('STA 30 scale');
set(gca,'xTick',[]);
set(gca,'yTick',[]);
%print(h,'-dpdf',sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/modelCellLNLN/model1/STA30.pdf'));

h=figure;
imagesc( repelem(sum(model.totalConeMap3D,3)==0,1,1,3).*ssf*20 + model.totalConeMap3D);
axis image
title('STA 20 scale');
set(gca,'xTick',[]);
set(gca,'yTick',[]);
%print(h,'-dpdf',sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/modelCellLNLN/model1/STA20.pdf'));
%% fit GMLM

 mask2 = logical(ones(size(movie,1),size(movie,2)));
 maskedMov= filterMov(movie,mask2,squeeze(model.ttf));
 maskedMov2=[maskedMov;ones(1,size(maskedMov,2))];
 binnedResponses = response';

 testLL_mc=[];
 dataLen_list = [500,750,1000,5000,10000,50000];
 for imc=1:20
     imc
 Tlen = size(maskedMov,2);
 times = randperm(Tlen);
 testtimes = times(1:5000);
 testLL_log =[];
 
 for itrainLen =  dataLen_list
     itrainLen
 traintimes = times(5001:5001+itrainLen);
 mm_train = maskedMov(:,traintimes); mm_test = maskedMov(:,testtimes);
 binnedResponses_train = binnedResponses(traintimes); binnedResponses_test = binnedResponses(testtimes);
 
 filteredStimDim =size(maskedMov,1);
 
 %  EM like Max Expected Likelihood .. 
 interval=1;
 fitGMLM1=cell(15,1); 
 nSU = 4;%1:15
 
 close all
 [fitGMLM,f_val] = fitGMLM_MEL_EM_bias(binnedResponses_train,mm_train,filteredStimDim,nSU,interval); 
 nTrails=1;
 [predictedResponse,lam,kx] = predictGMLM_bias_lr(fitGMLM,mm_test,nTrials,interval);lam=lam*120;
 testLL = (sum(lam)/120 - binnedResponses_test'*log(lam))/length(lam);
 testLL_log = [testLL_log;testLL];
 end
 
 testLL_mc(:,imc) =testLL_log;
 end

 
 meanLL = mean(testLL_mc,2);
 varLL =sqrt(var(testLL_mc,0,2)/size(testLL_mc,1));
 
 figure;
 ax = errorbar(log(dataLen_list),meanLL,varLL);
 