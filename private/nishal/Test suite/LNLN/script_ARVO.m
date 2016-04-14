
%% Make model 
model = model_LNLN();

%% Generate test response

pixelSz=8;
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

testLL_mc=[];
for imc=1:1
    imc
model = model_LNLN();
pixelSz=2;
model = setup_cutting_metric(model,pixelSz);

close all
% generate stimuli

sz = model.gridSzX/ (pixelSz*3);

Tlen = 120*60*180;
movie = (randn(sz,sz,Tlen)>0)-0.5;
dt=1/120;
nTrials=1;
[response,~] = generateResp_LNLN(model,movie,dt,nTrials);


 mask2 = logical(ones(size(movie,1),size(movie,2)));
 maskedMov= filterMov(movie,mask2,squeeze(model.ttf));
 maskedMov2=[maskedMov;ones(1,size(maskedMov,2))];
 binnedResponses = response';

 dataLen_list = [500,750,1000,5000,10000,50000,100000];

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

 
 meanLL = meduab(testLL_mc,2);
 varLL =sqrt(var(testLL_mc,0,2)/size(testLL_mc,1));
 
 figure;
 ax = errorbar(log(dataLen_list),meanLL,varLL);set(gca,'ylim',[-1 8]);
 
 
 %% 'cutting' metric
 
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

% fit sub-units

 mask2 = logical(ones(size(movie,1),size(movie,2)));
 maskedMov= filterMov(movie,mask2,squeeze(model.ttf));
 maskedMov2=[maskedMov;ones(1,size(maskedMov,2))];
 binnedResponses = response';
 Tlen = size(maskedMov,2);
 times = randperm(Tlen);
 testtimes = times(1:5000);
 traintimes = times(5001:end);
 mm_train = maskedMov(:,traintimes); mm_test = maskedMov(:,testtimes);
 binnedResponses_train = binnedResponses(traintimes); binnedResponses_test = binnedResponses(testtimes);
 filteredStimDim =size(maskedMov,1);
 interval=1;

 nSU=4%1:8
[fitGMLM,f_val] = fitGMLM_MEL_EM_bias(binnedResponses_train,mm_train,filteredStimDim,nSU,interval); 
 
sta_dim1 = size(mask2,1);
sta_dim2 = size(mask2,2);
indexedframe = reshape(1:sta_dim1*sta_dim2,[sta_dim1,sta_dim2]);
masked_frame = indexedframe(logical(mask2));

h=figure;
for ifilt=1:nSU
%subplot(szz(nSU,1),szz(nSU,2),ifilt)
subplot(1,nSU,ifilt);
u_spatial = reshape_vector(fitGMLM.Linear.filter{ifilt}(1:length(masked_frame)),masked_frame,indexedframe);

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

% plot weights
figure;
plot(model.SU_gang_weights.*sqrt(sum(model.cone_to_SU_connection,2)),'*')
title('SU to ganglion Weight distribution');

fitSUwt = [];
for isu=1:nSU
    fitSUwt = [fitSUwt ; norm(fitGMLM.Linear.filter{isu})];
end

figure;
plot(sort(fitSUwt,'descend'),'*');
title('Extracted SU magnitudes');


% radar chart 
dot_filter_su=zeros(nSU,model.nSU);
for ifilter = 1:nSU
    u_spatial = reshape_vector(fitGMLM.Linear.filter{ifilter}(1:length(masked_frame)),masked_frame,indexedframe);
    strongestFrame = u_spatial;
    szstr = size(strongestFrame,1);
    ssf = repelem(strongestFrame,model.gridSzX/szstr,model.gridSzX/szstr,1);
   [SU_inp]=filter_su_dot_LNLN(model,ssf);
    
    dot_filter_su(ifilter,:) = SU_inp'./sqrt(sum(model.cone_to_SU_connection.^2,2))'; % divide by sqrt(#cones) ? 
end

h =spider(dot_filter_su',sprintf('Extracted filters: %d',nSU));

% plot the messed up(permuted) sub-units

u_spatial_log = zeros(filteredStimDim,nSU);
h=figure;
for ifilt=1:nSU
%subplot(szz(nSU,1),szz(nSU,2),ifilt)

u_spatial_log(:,ifilt) =fitGMLM.Linear.filter{ifilt}(1:length(masked_frame));
end

for idim = 1:filteredStimDim
    perm = randperm(nSU);
u_spatial_log(idim,:)  = u_spatial_log(idim,perm);
end

for ifilt=1:nSU
    u_spatial = u_spatial_log(:,ifilt);
     u_spatial = reshape_vector(u_spatial,masked_frame,indexedframe);
subplot(1,nSU,ifilt);
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

%% compute parameters for the metric

%% compute value of cutting metric

for imc=1
    imc
    true_nSU = 12;
model = model_LNLN_parameterized(true_nSU);
pixelSz=16;
model = setup_cutting_metric(model,pixelSz);

close all
% generate stimuli

sz = model.gridSzX/ (pixelSz*3);

Tlen = 120*30*30;
movie = (randn(sz,sz,Tlen)>0)-0.5;
dt=1/120;
nTrials=1;
[response,~] = generateResp_LNLN(model,movie,dt,nTrials);


% fit sub-units
mask2 = logical(ones(size(movie,1),size(movie,2)));
sta_dim1 = size(mask2,1);
sta_dim2 = size(mask2,2);
indexedframe = reshape(1:sta_dim1*sta_dim2,[sta_dim1,sta_dim2]);
masked_frame = indexedframe(logical(mask2));

 maskedMov= filterMov(movie,mask2,squeeze(model.ttf));
 maskedMov2=[maskedMov;ones(1,size(maskedMov,2))];
 binnedResponses = response';

 Tlen = size(maskedMov,2);
 times = randperm(Tlen);
 testtimes = times(1:5000);
 traintimes = times(5001:end);
 mm_train = maskedMov(:,traintimes); mm_test = maskedMov(:,testtimes);
 binnedResponses_train = binnedResponses(traintimes); binnedResponses_test = binnedResponses(testtimes);
 filteredStimDim =size(maskedMov,1);
 
 %  EM like Max Expected Likelihood .. 
 interval=1;
 nSU = 4;
 close all
 [fitGMLM,f_val] = fitGMLM_MEL_EM_bias(binnedResponses_train,mm_train,filteredStimDim,nSU,interval); 
 nTrails=1;
 [predictedResponse,lam,kx] = predictGMLM_bias_lr(fitGMLM,mm_test,nTrials,interval);lam=lam*120;
 testLL = (sum(lam)/120 - binnedResponses_test'*log(lam))/length(lam);
 
u_spatial_log = zeros(filteredStimDim,nSU);

%h=figure;
for ifilt=1:nSU
%subplot(szz(nSU,1),szz(nSU,2),ifilt)
u_spatial_log(:,ifilt) =fitGMLM.Linear.filter{ifilt}(1:length(masked_frame));
end

u_spatial_log_perm = 0*u_spatial_log;
for idim = 1:filteredStimDim
    perm = randperm(nSU);
u_spatial_log_perm(idim,:)  = u_spatial_log(idim,perm);
end

[metric,metric_sus] = cutting_metric(model,u_spatial_log);

[metric_perm,metric_sus_perm] = cutting_metric(model,u_spatial_log_perm);


mc_data(imc).metric=metric;
mc_data(imc).metric_perm=metric_perm;
mc_data(imc).metric_sus=metric_sus;
mc_data(imc).metric_sus_perm=metric_sus_perm;
end

%%
m_log=[];m_perm_log=[];
for imc =1:length(mc_data)
m_log = [m_log;mc_data(imc).metric(1)];
m_perm_log = [m_perm_log;mc_data(imc).metric_perm(1)];
end


figure;
histogram(m_log,20);
hold on;%xlim([0,1]);

histogram(m_perm_log,20);
legend('fitting output','fitting output permuted');


m_log=[];m_perm_log=[];
for imc =1:length(mc_data)
m_log = [m_log;mc_data(imc).metric(2)];
m_perm_log = [m_perm_log;mc_data(imc).metric_perm(2)];
end


figure;
histogram(m_log,20);
hold on;%xlim([0,1]);

histogram(m_perm_log,20);
legend('fitting output','fitting output permuted');

%% plot u_spatial_log and u_spatial_log_perm

toplot = u_spatial_log;
for ifilt=1:nSU
    u_spatial = toplot(:,ifilt);
     u_spatial = reshape_vector(u_spatial,masked_frame,indexedframe);
subplot(1,nSU,ifilt);
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


[metric,metric_sus] = cutting_metric(model,toplot);
suptitle(sprintf('Metrics 1: %0.02f, 2: %0.02f',metric(1),metric(2)));

%% calculate cutting metric baseline
%noise=21;
metric_log=[];
for imc=1:100
imc
    true_nSU = 12;
model = model_LNLN_parameterized(true_nSU);

%model = model_LNLN();
pixelSz=16;
model = setup_cutting_metric(model,pixelSz);

close all
% generate stimuli
% 
% sz = model.gridSzX/ (pixelSz*3);
% 
% Tlen = 120*30*30;
% movie = (randn(sz,sz,Tlen)>0)-0.5;
% dt=1/120;
% nTrials=1;
% [response,~] = generateResp_LNLN(model,movie,dt,nTrials);
% 
% % calculate STA
% idx = 1:Tlen;
% spktm = idx(response~=0 & idx >30);
% STA= zeros(sz,sz,30);
% for itime=spktm
% STA = STA + movie(:,:,itime-29:itime);
% end
% STA = STA / length(spktm);
% 
% figure;
% for itime=26
%     itime
% imagesc(STA(:,:,itime));
% colormap gray
% caxis([min(STA(:)),max(STA(:))]);
% %pause
% end
% 
% pixx = STA(end-3:end,end-3:end,:);
% noise = sqrt(var(pixx(:)));

% make 'good sub-units'
% select number of true su in 'good' fitted SU
nSU_fit = 4;
nSU_true = model.nSU;
ntrue_per_fit = ceil(nSU_true/nSU_fit);
togo=1;
while togo==1
    nSUi=[];
    for isu =[1:nSU_fit-1]
        nSUi = [nSUi;poissrnd(ntrue_per_fit-1)];
    end
    nSUi = [nSUi;nSU_true-sum(nSUi)];
    nSUi(nSUi<=0)=1;
    if(sum(nSUi)==true_nSU)
        togo=0;
    end
end

% evenly spread the SU
%nSUi = ones(nSU_fit,1)*ntrue_per_fit;

% randomly (in future, prefer near ones?) select sub-units
select_matrix =zeros(nSU_fit,nSU_true);
for isu=1:nSU_fit
select_matrix(isu,:) = [zeros(1,sum(nSUi(1:isu-1))),ones(1,nSUi(isu)),zeros(1,nSU_true-sum(nSUi(1:isu)))];
end
perms = randperm(nSU_true);
select_matrix = select_matrix(:,perms);

u_spatial_log = (select_matrix * model.su_lowres')';
noise = norm(sum(u_spatial_log,2))/(10*sqrt(sum(abs(sum(u_spatial_log,2))>0.01)));

u_spatial_log = u_spatial_log+ noise*randn(size(u_spatial_log));
[metric,metric_sus] = cutting_metric(model,u_spatial_log);
nSU = nSU_fit

metric_log = [metric_log;metric];

end

figure;
histogram(metric_log(:,1))

%% plot histograms
data = load('~/Google Drive/Presentations/Presentations/ASM_figures/modelCell_metric.mat');
baseline = load('~/Google Drive/Presentations/Presentations/ASM_figures/modelCell_metric_baseline2.mat');
mc_data = data.mc_data;
metric_log = baseline.metric_log;

m_log=[];m_perm_log=[];
for imc =1:length(mc_data)
m_log = [m_log;mc_data(imc).metric(1)];
m_perm_log = [m_perm_log;mc_data(imc).metric_perm(1)];
end


figure;
histogram(m_log,20);
hold on;%xlim([0,1]);

histogram(m_perm_log,20);
hold on;
histogram(metric_log(:,1),20);
legend('fitting output','fitting output permuted','baseline');


m_log=[];m_perm_log=[];
for imc =1:length(mc_data)
m_log = [m_log;mc_data(imc).metric(2)];
m_perm_log = [m_perm_log;mc_data(imc).metric_perm(2)];
end


figure;
histogram(m_log,20);
hold on;%xlim([0,1]);

histogram(m_perm_log,20);

hold on;
histogram(metric_log(:,2),20);
legend('fitting output','fitting output permuted','baseline');


%% fine Res , long WN - LNL spikes, fit
% as its fine resolution, we do not need to take into account stixels! -
% excite cones independently with gaussian noise.

for imc=1:100
    imc
model = model_LNLN_fineRes();

% Generate test response


Tlen = 120*10;
conemovie = 25*(randn(Tlen,model.nCones));
dt=1/120;
nTrials=30;
[response,su_inp,conefilteredMov] = generateResp_LNLN_coneMov(model,conemovie,dt,nTrials);
h= figure;
plotSpikeRaster(response~=0,'PlotType','vertline');
%print(h,'-dpdf',sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/modelCellLNLN/model1/sample_firing.pdf'));

Tlen = 120*60*180;
conemovie = 25*(randn(Tlen,model.nCones));
dt=1/120;
nTrials=1;
[response,~,conefilteredMov] = generateResp_LNLN_coneMov(model,conemovie,dt,nTrials);

% cone sub-units. 

maskedMov = conefilteredMov/1000;
binnedResponses = response';
filteredStimDim =size(maskedMov,1);
interval=1; 
nSU=12;
[fitGMLM,f_val] = fitGMLM_MEL_EM_bias(binnedResponses,maskedMov,filteredStimDim,nSU,interval); 



u_spatial_log = zeros(filteredStimDim,nSU);

%h=figure;
for ifilt=1:nSU
%subplot(szz(nSU,1),szz(nSU,2),ifilt)
u_spatial_log(:,ifilt) =fitGMLM.Linear.filter{ifilt}(1:length(masked_frame));
end

u_spatial_log_perm = 0*u_spatial_log;
for idim = 1:filteredStimDim
    perm = randperm(nSU);
u_spatial_log_perm(idim,:)  = u_spatial_log(idim,perm);
end

[metric,metric_sus] = cutting_metric(model,u_spatial_log);

[metric_perm,metric_sus_perm] = cutting_metric(model,u_spatial_log_perm);


mc_data(imc).metric=metric;
mc_data(imc).metric_perm=metric_perm;
mc_data(imc).metric_sus=metric_sus;
mc_data(imc).metric_sus_perm=metric_sus_perm;
end



% plot sub-units
cols = distinguishable_colors(nSU+1);
figure;
for isu_fitted=1:nSU
    subplot(3,4,isu_fitted);
    for isu_true = 1:model.nSU
        cones = (model.cone_su_idx==isu_true);
        scatter(model.conesX(cones),model.conesY(cones),7*abs(fitGMLM.Linear.filter{isu_fitted}(cones)),cols(isu_true,:),'filled');
        hold on;
        xlim([min(model.conesX),max(model.conesX)]);
        ylim([min(model.conesY),max(model.conesY)]);
        axis square
        set(gca,'xTick',[]);set(gca,'yTick',[]);
        title(sprintf('fitted SU#',isu_fitted));
    end
     set(gca,'visible','off');
end
%

figure;
for isu_true = 1:model.nSU
        cones = (model.cone_su_idx==isu_true);
        scatter(model.conesX(cones),model.conesY(cones),250,cols(isu_true,:),'filled');
        hold on;
        xlim([min(model.conesX),max(model.conesX)]);
        ylim([min(model.conesY),max(model.conesY)]);
        axis square
        set(gca,'xTick',[]);set(gca,'yTick',[]);
end
 set(gca,'visible','off');
 
 %% fine Res , long WN - LNLN spikes, LNL fit
% as its fine resolution, we do not need to take into account stixels! -
% excite cones independently with gaussian noise.

model = model_LNLN_fineRes_fancyNL();

% Generate test response
Tlen = 120*10;
conemovie = 25*(randn(Tlen,model.nCones));
dt=1/120;
nTrials=30;
[response,su_inp,conefilteredMov] = generateResp_LNLN_coneMov_nl(model,conemovie,dt,nTrials);
h= figure;
plotSpikeRaster(response~=0,'PlotType','vertline');
%print(h,'-dpdf',sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/modelCellLNLN/model1/sample_firing.pdf'));

Tlen = 120*60*30;
conemovie = 25*(randn(Tlen,model.nCones));
dt=1/120;
nTrials=1;
[response,~,conefilteredMov] = generateResp_LNLN_coneMov_nl(model,conemovie,dt,nTrials);

% cone sub-units. 

maskedMov = conefilteredMov/1000;
binnedResponses = response';
filteredStimDim =size(maskedMov,1);
interval=1; 
nSU=12;
[fitGMLM,f_val] = fitGMLM_MEL_EM_bias(binnedResponses,maskedMov,filteredStimDim,nSU,interval); 

% plot sub-units
cols = distinguishable_colors(nSU+1);
figure;
for isu_fitted=1:nSU
    subplot(3,4,isu_fitted);
    for isu_true = 1:model.nSU
        cones = (model.cone_su_idx==isu_true);
        scatter(model.conesX(cones),model.conesY(cones),7*abs(fitGMLM.Linear.filter{isu_fitted}(cones)),cols(isu_true,:),'filled');
        hold on;
        xlim([min(model.conesX),max(model.conesX)]);
        ylim([min(model.conesY),max(model.conesY)]);
        axis square
        set(gca,'xTick',[]);set(gca,'yTick',[]);
        title(sprintf('fitted SU#',isu_fitted));
    end
     set(gca,'visible','off');
end
%

figure;
for isu_true = 1:model.nSU
        cones = (model.cone_su_idx==isu_true);
        scatter(model.conesX(cones),model.conesY(cones),250,cols(isu_true,:),'filled');
        hold on;
        xlim([min(model.conesX),max(model.conesX)]);
        ylim([min(model.conesY),max(model.conesY)]);
        axis square
        set(gca,'xTick',[]);set(gca,'yTick',[]);
end
 set(gca,'visible','off');
 
 
 %% smoothness/spatial locality metric
 
for imc=1:10
    imc
    true_nSU = 12;
model = model_LNLN_parameterized(true_nSU);
pixelSz=16;
model = setup_cutting_metric(model,pixelSz);

close all
% generate stimuli

sz = model.gridSzX/ (pixelSz*3);

Tlen = 120*30*30;
movie = (randn(sz,sz,Tlen)>0)-0.5;
dt=1/120;
nTrials=1;
[response,~] = generateResp_LNLN(model,movie,dt,nTrials);


% fit sub-units
mask2 = logical(ones(size(movie,1),size(movie,2)));
sta_dim1 = size(mask2,1);
sta_dim2 = size(mask2,2);
indexedframe = reshape(1:sta_dim1*sta_dim2,[sta_dim1,sta_dim2]);
masked_frame = indexedframe(logical(mask2));

 maskedMov= filterMov(movie,mask2,squeeze(model.ttf));
 maskedMov2=[maskedMov;ones(1,size(maskedMov,2))];
 binnedResponses = response';

 Tlen = size(maskedMov,2);
 times = randperm(Tlen);
 testtimes = times(1:5000);
 traintimes = times(5001:end);
 mm_train = maskedMov(:,traintimes); mm_test = maskedMov(:,testtimes);
 binnedResponses_train = binnedResponses(traintimes); binnedResponses_test = binnedResponses(testtimes);
 filteredStimDim =size(maskedMov,1);
 
 %  EM like Max Expected Likelihood .. 
 interval=1;
 nSU = 4;
 close all
 [fitGMLM,f_val] = fitGMLM_MEL_EM_bias(binnedResponses_train,mm_train,filteredStimDim,nSU,interval); 
 nTrails=1;
 [predictedResponse,lam,kx] = predictGMLM_bias_lr(fitGMLM,mm_test,nTrials,interval);lam=lam*120;
 testLL = (sum(lam)/120 - binnedResponses_test'*log(lam))/length(lam);
 
u_spatial_log = zeros(filteredStimDim,nSU);

%h=figure;
for ifilt=1:nSU
%subplot(szz(nSU,1),szz(nSU,2),ifilt)
u_spatial_log(:,ifilt) =fitGMLM.Linear.filter{ifilt}(1:length(masked_frame));
end

u_spatial_log_perm = 0*u_spatial_log;
for idim = 1:filteredStimDim
    perm = randperm(nSU);
u_spatial_log_perm(idim,:)  = u_spatial_log(idim,perm);
end

spatial_scale=2;
[metric,metric_sus] = smoothness_metric(u_spatial_log,logical(ones(10,10)),spatial_scale);
[metric_perm,metric_sus_perm] = smoothness_metric(u_spatial_log_perm,logical(ones(10,10)),spatial_scale);

[metric,metric_sus] = smoothness_metric2(u_spatial_log,logical(ones(10,10)),spatial_scale);
[metric_perm,metric_sus_perm] = smoothness_metric2(u_spatial_log_perm,logical(ones(10,10)),spatial_scale);


mc_data(imc).metric=metric;
mc_data(imc).metric_perm=metric_perm;
mc_data(imc).metric_sus=metric_sus;
mc_data(imc).metric_sus_perm=metric_sus_perm;
end

metric=[];metric_perm=[];
for imc=1:length(mc_data)
metric = [metric;mc_data(imc).metric];
metric_perm = [metric_perm;mc_data(imc).metric_perm];
end
 
figure;
histogram(metric,10);
hold on;
histogram(metric_perm,10);