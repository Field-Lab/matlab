
%% Make model 
model = model_LNLN();

%% Generate test response
sz = model.gridSzX/ (16*3);


Tlen = 120*10;
movie = (randn(sz,sz,Tlen)>0)-0.5;
dt=1/120;
nTrials=30;
[response,~] = generateResp_LNLN(model,movie,dt,nTrials);
h= figure;
plotSpikeRaster(response~=0,'PlotType','vertline');
%print(h,'-dpdf',sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/modelCellLNLN/model1/sample_firing.pdf'));

Tlen = 120*60*30;
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
for itime=1:30
    itime
imagesc(STA(:,:,itime));
colormap gray
caxis([min(STA(:)),max(STA(:))]);
pause
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
print(h,'-dpdf',sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/modelCellLNLN/model1/STA30.pdf'));

h=figure;
imagesc( repelem(sum(model.totalConeMap3D,3)==0,1,1,3).*ssf*20 + model.totalConeMap3D);
axis image
title('STA 20 scale');
set(gca,'xTick',[]);
set(gca,'yTick',[]);
print(h,'-dpdf',sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/modelCellLNLN/model1/STA20.pdf'));
%% fit GMLM

 mask2 = logical(ones(size(movie,1),size(movie,2)));
 maskedMov= filterMov(movie,mask2,squeeze(model.ttf));
 maskedMov2=[maskedMov;ones(1,size(maskedMov,2))];

 binnedResponses = response';

 filteredStimDim =size(maskedMov,1);
 
 %  EM like Max Expected Likelihood .. 
 interval=1;
 %[fitGMLM,output] = fitGMLM_MEL_EM(binnedResponses,maskedMov2,8,4,interval);   
 
 fitGMLM1=cell(15,1); 
 for nSU=1:15
 nSU
 fitGMLM_log = cell(50,1);
 fval_log = zeros(50,1);
 for ifit = 1:50
     ifit
     close all
 [fitGMLM,f_val] = fitGMLM_MEL_EM_bias(binnedResponses,maskedMov,filteredStimDim,nSU,interval); 
 
 fitGMLM_log{ifit} = fitGMLM;
 fval_log(ifit) = f_val;
 
 end
 fitGMLM1{nSU} = fitGMLM_log;
 save(sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/modelCellLNLN/model1/stix 16, 90 min/fit_nSU_%d.mat',nSU),'fitGMLM_log','fval_log','model','mask2');
 end
  

%% Show filters
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
    4,4];

for nSU=4%1:8
 fitGMLM{1} =fitGMLM1{7};   
sta_dim1 = size(mask2,1);
sta_dim2 = size(mask2,2);
indexedframe = reshape(1:sta_dim1*sta_dim2,[sta_dim1,sta_dim2]);
masked_frame = indexedframe(logical(mask2));

h=figure;
for ifilt=1:nSU
%subplot(szz(nSU,1),szz(nSU,2),ifilt)
subplot(1,nSU,ifilt);
u_spatial = reshape_vector(fitGMLM{1}.Linear.filter{ifilt}(1:length(masked_frame)),masked_frame,indexedframe);

strongestFrame = -1*u_spatial/max(abs(u_spatial(:)))+0.5;
szstr = size(strongestFrame,1);
ssf = repelem(strongestFrame,model.gridSzX/szstr,model.gridSzX/szstr,3);
mmask = sum(model.totalConeMap3D,3)==0;
xx = repelem(mmask,1,1,3).*ssf; 
aa = model.totalConeMap3D;  aa(repelem(sum(aa,3),1,1,3)<0.5)=aa(repelem(sum(aa,3)<0.5,1,1,3))+0.5;aa =(aa-0.5)/max(aa(:)-0.5) + 0.5;
mxt = ((aa(repelem(sum(aa,3)~=1.5,1,1,3))));mx= max(mxt(:));
xx(repelem(sum(aa,3)~=1.5,1,1,3))= mxt;
xx = xx(1:ceil(max(model.conesX))+40,1:ceil(max(model.conesY))+40,:);

% mag = repelem(sum(model.totalConeMap3D,3)==0,1,1,3).*ssf + model.totalConeMap3D*2;
% mag = mag(1:ceil(max(model.conesX))+40,1:ceil(max(model.conesY))+40,:);

imagesc(xx);
axis image
set(gca,'xTick',[]);
set(gca,'yTick',[]);

title(sprintf('SU # %d',ifilt));
%caxis([-max(mag(:)),max(mag(:))]);
%caxis([0,max(aa(:))])
end

%print(h,'-dpdf',sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/modelCellLNLN/model1/Exp_ASM_MEL_EM_filters SU%d.pdf',nSU));

% plot weights
figure;
plot(sort(model.SU_gang_weights.*sqrt(sum(model.cone_to_SU_connection,2)),'descend'),'*')
title('SU to ganglion Weight distribution');

fitSUwt = [];
for isu=1:nSU
    fitSUwt = [fitSUwt ; norm(fitGMLM{1}.Linear.filter{isu})];
end

figure;
plot(sort(fitSUwt,'descend'),'*');
title('Extracted SU magnitudes');


% radar chart 
dot_filter_su=zeros(nSU,model.nSU);
for ifilter = 1:nSU
    u_spatial = reshape_vector(fitGMLM{1}.Linear.filter{ifilter}(1:length(masked_frame)),masked_frame,indexedframe);
    strongestFrame = u_spatial;
    szstr = size(strongestFrame,1);
    ssf = repelem(strongestFrame,model.gridSzX/szstr,model.gridSzX/szstr,1);
   [SU_inp]=filter_su_dot_LNLN(model,ssf);
    
    dot_filter_su(ifilter,:) = SU_inp'./sqrt(sum(model.cone_to_SU_connection.^2,2))'; % divide by sqrt(#cones) ? 
end

h =spider(dot_filter_su',sprintf('Extracted filters: %d',nSU));

%print(h,'-dpdf',sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/modelCellLNLN/model1/spider SU%d.pdf',nSU));
end

%% Weight magnitudes

ratioFilters=[];
maxMag = [];
for nSU = 1:8
fitGMLM  = fitGMLM1{nSU};

filterMag = [];
for ifilter=1:nSU
filterMag = [filterMag;norm(fitGMLM{1}.Linear.filter{ifilter})];
end

maxMag(nSU) = max(filterMag);
ratioFilters(nSU) = min(filterMag)/max(filterMag);

end

figure;
plot(1:8,ratioFilters,'-*');
ylim([0,1]);
title('Min/Max filter magnitude v/s nSU extracted');


figure;
plot(1:8,maxMag,'-*');
%ylim([0,1]);
title('Max filter magnitude v/s nSU extracted');

%% STA null stimulus .. 


movieLen=120*12;
stas{1}=zeros(size(STA,1),size(STA,2),3,size(STA,3));
stas{1}(:,:,1,:)=STA;
stas{1}(:,:,2,:)=STA;
stas{1}(:,:,3,:)=STA;
cell_params.STAlen=14;
[new_stas,mask,CellMasks]=clipSTAs(stas,cell_params)
Filtdim1=size(STA,1);Filtdim2=size(STA,2);Filtlen=30;
null_compute_subUnit_test
 
nTrials=50;
dt=1/120;
[respOrig,~] = generateResp_LNLN(model,mov_orig2,dt,nTrials);
[respNull,~] = generateResp_LNLN(model,mov_new2,dt,nTrials);
[x1,y1]=plotSpikeRaster(logical(respOrig));
[x2,y2]=plotSpikeRaster(logical(respNull));

figure;
plot(x1,y1+nTrials,'k');
hold on;
plot(x2,y2,'r');
xlim([0,max(x2(:))]);
%% Predict null response

R2_orig=zeros(15,1);
R2_null=zeros(15,1);

movie2=mov_orig2;
mask2 = logical(ones(size(movie2,1),size(movie2,2)));
maskedMov_orig= filterMov(movie2,mask2,squeeze(model.ttf));

% Predict response for null movie
movie2=mov_new2; 
mask2 = logical(ones(size(movie2,1),size(movie2,2)));
maskedMov_null= filterMov(movie2,mask2,squeeze(model.ttf));


for nSU=1:15
    
    load(sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/modelCellLNLN/model1/stix 16, 90 min/fit_nSU_%d.mat',nSU))
    
    for ifit =3%1:50
    fitGMLM = fitGMLM_log{ifit};
% Predict response for original movie

filteredStimDim =size(maskedMov_orig,1);
nTrials=50;
interval=1;
[predictedRespOrig,lamOrig] = predictGMLM_bias_lr(fitGMLM,maskedMov_orig,nTrials,interval);
predictedRespOrig=predictedRespOrig';

%
filteredStimDim =size(maskedMov_null,1);
nTrials=50;
interval=1;
[predictedRespNull,lamNull] = predictGMLM_bias_lr(fitGMLM,maskedMov_null,nTrials,interval);
predictedRespNull=predictedRespNull';

[x1,y1] = plotSpikeRaster(logical(respOrig));
[x2,y2] = plotSpikeRaster(predictedRespOrig>0);
[x3,y3] = plotSpikeRaster(logical(respNull));
[x4,y4] = plotSpikeRaster(predictedRespNull>0);
% 
figure;
plot(x1,y1+50*3,'k');
hold on;
plot(x2,y2+50*2,'r');
hold on;
plot(x3,y3+50*1,'k');
hold on;
plot(x4,y4,'r');


% % Calculate R-2 value
%    convolve=150;
%   
%     binSz=1/120;
%     R2_log=[];
%     
%     % Original real
%     len = size(respOrig,2);
%     [PSTHOrig,time]=calculate_psth_fcn2(convolve,binSz,len,respOrig);
%     
%     % Original predicted
%     len = size(predictedRespOrig,2);
%     [PSTHpredOrig,time]=calculate_psth_fcn2(convolve,binSz,len,predictedRespOrig);
%      
%       % Null real
%     len = size(respNull,2);
%     [PSTHNull,time]=calculate_psth_fcn2(convolve,binSz,len,respNull);
%     
%     % Null predicted
%     len = size(predictedRespNull,2);
%     [PSTHpredNull,time]=calculate_psth_fcn2(convolve,binSz,len,predictedRespNull);
    
%     figure;
%     plot(time,PSTHOrig);hold on;
%     plot(time,PSTHpredOrig);hold on;
%     plot(time,PSTHNull);hold on;
%     plot(time,PSTHpredNull);hold on;
%     
%     idx = time>0.3*max(time) & time<0.7*max(time); % remove starting and ending.
%    R2_orig(nSU) = R_2_value(PSTHOrig(idx)',PSTHpredOrig(idx)') 
%    R2_null(nSU) = R_2_value(PSTHNull(idx)',PSTHpredNull(idx)') 
    end
end

h=figure;
plot(R2_orig,'--*');
hold on
plot(R2_null,'--*');
xlabel('number of S.U.');
ylabel('R2 value');
ylim([0,1]);
xlim([0,15]);
title('Prediction accuracy v/s number of SU')
legend('Original movie R2','Null movie R2');
print(h,'-dpdf',sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/modelCellLNLN/model1/Rsqvalue_vs_SU.pdf'));


%%
%predict null response

R2_log_orig=zeros(15,50);
R2_log_null=zeros(15,50);

movie2=mov_orig2;
mask2 = logical(ones(size(movie2,1),size(movie2,2)));
maskedMov_orig= filterMov(movie2,mask2,squeeze(model.ttf));

% Predict response for null movie
movie2=mov_new2; 
mask2 = logical(ones(size(movie2,1),size(movie2,2)));
maskedMov_null= filterMov(movie2,mask2,squeeze(model.ttf));


for nSU=1:8
    nSU
    load(sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/modelCellLNLN/model1/stix 16, 90 min/fit_nSU_%d.mat',nSU))
    
    for ifit =1:50
        ifit
        fitGMLM = fitGMLM_log{ifit};
% Predict response for original movie

filteredStimDim =size(maskedMov_orig,1);
nTrials=50;
interval=1;
[predictedRespOrig,lamOrig] = predictGMLM_bias_lr(fitGMLM,maskedMov_orig,nTrials,interval);
predictedRespOrig=predictedRespOrig';

%
filteredStimDim =size(maskedMov_null,1);
nTrials=50;
interval=1;
[predictedRespNull,lamNull] = predictGMLM_bias_lr(fitGMLM,maskedMov_null,nTrials,interval);
predictedRespNull=predictedRespNull';


v = R_2_value(lamOrig,respOrig(1,:)');
R2_log_orig(nSU,ifit) =v;

v = R_2_value(lamNull,respNull(1,:)');
R2_log_null(nSU,ifit) =v;
    end
end

% plot stuff

r2_su = R2_log_orig;
r2_su2 = R2_log_null;
nc=15;
       h=figure;
       errorbar(1:nc,mean(r2_su,2),sqrt(var(r2_su')),'-*');
       hold on;
       plot(1:nc,max(r2_su'),'-*');
       hold on;
       errorbar(1:nc,mean(r2_su2,2),sqrt(var(r2_su2')),'-*');
       hold on;
       plot(1:nc,max(r2_su2'),'-*');
       hold on;
       
       legend('mean for BW movie','max for BW movie','mean for Null movie','max for Null movie','location','best');
       title('R^2 value v/s max. number of sub-units');
       xlabel('Number of sub-units');
       ylabel('R^2');
 hgexport(h,sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/modelCellLNLN/model1/stix 16, 90 min/R2_Null.eps'));
 
 %%
 
movie2=mov_orig2;
mask2 = logical(ones(size(movie2,1),size(movie2,2)));
maskedMov_orig= filterMov(movie2,mask2,squeeze(model.ttf));

% Predict response for null movie
movie2=mov_new2; 
mask2 = logical(ones(size(movie2,1),size(movie2,2)));
maskedMov_null= filterMov(movie2,mask2,squeeze(model.ttf));

nSU=4; ifit=3;
load(sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/modelCellLNLN/model1/stix 16, 90 min/fit_nSU_%d.mat',nSU))

fitGMLM = fitGMLM_log{ifit};

filteredStimDim =size(maskedMov_null,1);
nTrials=50;
interval=1;
[predictedRespNull,lamNull] = predictGMLM_bias_lr(fitGMLM,maskedMov_orig,nTrials,interval);
predictedRespNull=predictedRespNull';

[x3,y3] = plotSpikeRaster(logical(respNull));
[x4,y4] = plotSpikeRaster(predictedRespNull>0);

nSU=1;ifit=3;

load(sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/modelCellLNLN/model1/stix 16, 90 min/fit_nSU_%d.mat',nSU))

fitGMLM = fitGMLM_log{ifit};

filteredStimDim =size(maskedMov_null,1);
nTrials=50;
interval=1;
[predictedRespNull,lamNull] = predictGMLM_bias_lr(fitGMLM,maskedMov_orig,nTrials,interval);
predictedRespNull=predictedRespNull';

[x5,y5] = plotSpikeRaster(predictedRespNull>0);

figure;
plot(x4/120,y4+50*2,'r');
hold on;
plot(x3/120,y3+50*1,'k');
hold on;
plot(x5/120,y5+50*0,'r');

set(gca,'yTick',[]);
xlabel('Time(s)');
xlim([2,12]);


%% long null stimulus .. 

% generate response to long null stimulus

movieLen=120*60*30;
stas{1}=zeros(size(STA,1),size(STA,2),3,size(STA,3));
stas{1}(:,:,1,:)=STA;
stas{1}(:,:,2,:)=STA;
stas{1}(:,:,3,:)=STA;
cell_params.STAlen=14;
[new_stas,mask,CellMasks]=clipSTAs(stas,cell_params)
Filtdim1=size(STA,1);Filtdim2=size(STA,2);Filtlen=30;


% null_compute_subUnit_test
model_spatial_nulling


nTrials=1;
dt=1/120;
[respOrig,~] = generateResp_LNLN(model,mov_orig2,dt,nTrials);
[respNull,~] = generateResp_LNLN(model,mov_new2,dt,nTrials);

movie = mov_orig2; response = respOrig;
model_GMLM_fit;
fitGMLM_log_orig =  fitGMLM1;


nSU=3; 
mask2 = logical(ones(size(movie,1),size(movie,2)));
initialData.initialFilters = 2*(rand(sum(double(mask2(:)))*nSU,1)-0.5);
initialData.initalbias = 2*(rand(nSU,1)-0.5);

% analyze with different training length
movie = mov_orig2; response = respOrig;
[fitGMLM_log_ori,mse_data_ori] = model_GMLM_fit_diff_data_len(movie,response,model,nSU,initialData,mask2,P_mat);

% analyze with different training length
movie = mov_new2; response = respNull;
[fitGMLM_log_null,mse_data_null] = model_GMLM_fit_diff_data_len(movie,response,model,nSU,initialData,mask2,P_mat);




