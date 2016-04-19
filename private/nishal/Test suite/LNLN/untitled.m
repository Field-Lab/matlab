
  
model = model_LNLN();
pixelSz=16;
model = setup_cutting_metric(model,pixelSz);

close all
% generate stimuli

sz = model.gridSzX/ (pixelSz*3);
% % % get sig_stixels 
Tlen = 120*30*30;
movie = (randn(sz,sz,Tlen)>0)-0.5;
dt=1/120;
nTrials=1;
[response,~] = generateResp_LNLN(model,movie,dt,nTrials);

% Calulate STA 
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

% find significant stixels
STA3D = zeros(size(STA,1),size(STA,2),3,size(STA,3));
STA3D(:,:,1,:)=STA;
STA3D(:,:,2,:)=STA;
STA3D(:,:,3,:)=STA;
[sig_stixels, params, rf_strength_out] = significant_stixels(STA3D);



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
 
u_spatial_log_fitted = zeros(filteredStimDim,nSU);

%h=figure;
for ifilt=1:nSU
%subplot(szz(nSU,1),szz(nSU,2),ifilt)
u_spatial_log_fitted(:,ifilt) =fitGMLM.Linear.filter{ifilt}(1:length(masked_frame));
end



% gaussian partitioning!

[xgrid,ygrid] = meshgrid(1:sz,1:sz);

xx = xgrid(sig_stixels);
yy = ygrid(sig_stixels);
nSU=4;sig = sqrt((max(xx) - min(xx))*(max(yy) - min(yy)) / nSU)/4;

[idx,C] = kmeans([xx,yy],nSU);
filteredStimDim = numel(STA(:,:,1));
u_spatial_log = zeros(filteredStimDim,nSU);


% resample cluster centers and reassign cluster identities.
resampleC = [];
for isu=1:nSU
xxsu = xx(idx==isu);
yysu = yy(idx==isu);
per= randperm(numel(xxsu));
resampleC = [resampleC;xxsu(per(1)),yysu(per(1))];
end

dist = zeros(numel(xx),nSU);
for isu=1:nSU
    xy_center = resampleC(isu,:);
    dist(:,isu) = sqrt((xx - xy_center(1)).^2 + (yy-xy_center(2)).^2);
end
[v,idx] = min(dist,[],2);

% figure; 
% for isu=1:nSU
%     xxsu = xx(idx==isu); yysu = yy(idx==isu);
%     sus = zeros(sz,sz);
%     for ix=1:numel(xxsu)
%     sus(yysu(ix),xxsu(ix))=1;%normpdf(v(xxsu(ix)==xx & yysu(ix)==yy),0,sig);
%     end
%     u_spatial_log(:,isu) = sus(:);
%     subplot(1,nSU,isu);
% imagesc(sus);axis image
% end


figure; 
for isu=1:nSU
    xxsu = xx; yysu = yy;
    sus = zeros(sz,sz);
    for ix=1:numel(xxsu)
    sus(yysu(ix),xxsu(ix))=normpdf(dist(xxsu(ix)==xx & yysu(ix)==yy,isu),0,sig) ;
    end
    u_spatial_log(:,isu) = sus(:);
    subplot(1,nSU,isu);
imagesc(sus);axis image
end

noise = norm(sum(u_spatial_log,2))/(10*sqrt(sum(abs(sum(u_spatial_log,2))>0.01)));

u_spatial_log_gauss = u_spatial_log+ noise*randn(size(u_spatial_log));


%%
for imc=1:100
    imc
   
model = model_LNLN();
pixelSz=16;
model = setup_cutting_metric(model,pixelSz);

close all
% generate stimuli

sz = model.gridSzX/ (pixelSz*3);
% % % get sig_stixels 
% Tlen = 120*30*30;
% movie = (randn(sz,sz,Tlen)>0)-0.5;
% dt=1/120;
% nTrials=1;
% [response,~] = generateResp_LNLN(model,movie,dt,nTrials);
% 
% % Calulate STA 
% idx = 1:Tlen;
% spktm = idx(response~=0 & idx >30);
% 
% STA= zeros(sz,sz,30);
% 
% 
% for itime=spktm
% STA = STA + movie(:,:,itime-29:itime);
% end
% 
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
% % find significant stixels
% STA3D = zeros(size(STA,1),size(STA,2),3,size(STA,3));
% STA3D(:,:,1,:)=STA;
% STA3D(:,:,2,:)=STA;
% STA3D(:,:,3,:)=STA;
% [sig_stixels, params, rf_strength_out] = significant_stixels(STA3D);
[xgrid,ygrid] = meshgrid(1:sz,1:sz);

xx = xgrid(sig_stixels);
yy = ygrid(sig_stixels);
nSU=4;sig = sqrt((max(xx) - min(xx))*(max(yy) - min(yy)) / nSU)/2;

[idx,C] = kmeans([xx,yy],nSU);
filteredStimDim = numel(STA(:,:,1));
u_spatial_log = zeros(filteredStimDim,nSU);

% resample cluster centers and reassign cluster identities.
resampleC = [];
for isu=1:nSU
xxsu = xx(idx==isu);
yysu = yy(idx==isu);
per= randperm(numel(xxsu));
resampleC = [resampleC;xxsu(per(1)),yysu(per(1))];
end

dist = zeros(numel(xx),nSU);
for isu=1:nSU
    xy_center = resampleC(isu,:);
    dist(:,isu) = sqrt((xx - xy_center(1)).^2 + (yy-xy_center(2)).^2);
end
[v,idx] = min(dist,[],2);

% figure; 
% for isu=1:nSU
%     xxsu = xx(idx==isu); yysu = yy(idx==isu);
%     sus = zeros(sz,sz);
%     for ix=1:numel(xxsu)
%     sus(yysu(ix),xxsu(ix))=1;%normpdf(v(xxsu(ix)==xx & yysu(ix)==yy),0,sig);
%     end
%     u_spatial_log(:,isu) = sus(:);
%     subplot(1,nSU,isu);
% imagesc(sus);axis image
% end


figure; 
for isu=1:nSU
    xxsu = xx; yysu = yy;
    sus = zeros(sz,sz);
    for ix=1:numel(xxsu)
    sus(yysu(ix),xxsu(ix))=normpdf(dist(xxsu(ix)==xx & yysu(ix)==yy,isu),0,sig) ;
    end
    u_spatial_log(:,isu) = sus(:);
    subplot(1,nSU,isu);
imagesc(sus);axis image
end

noise = norm(sum(u_spatial_log,2))/(10*sqrt(sum(abs(sum(u_spatial_log,2))>0.01)));

u_spatial_log = u_spatial_log+ noise*randn(size(u_spatial_log));

[metric,metric_sus] = cutting_metric(model,u_spatial_log);

%[metric_perm,metric_sus_perm] = cutting_metric(model,u_spatial_log_perm);


mc_data(imc).metric=metric;
%mc_data(imc).metric_perm=metric_perm;
mc_data(imc).metric_sus=metric_sus;
%mc_data(imc).metric_sus_perm=metric_sus_perm;
end

%%
m_log=[];m_perm_log=[];
for imc =1:length(mc_data)
m_log = [m_log;mc_data(imc).metric(1)];
%m_perm_log = [m_perm_log;mc_data(imc).metric_perm(1)];
end


figure;
histogram(m_log,20);
hold on;%xlim([0,1]);

%histogram(m_perm_log,20);
legend('fitting output');


m_log=[];m_perm_log=[];
for imc =1:length(mc_data)
m_log = [m_log;mc_data(imc).metric(2)];
%m_perm_log = [m_perm_log;mc_data(imc).metric_perm(2)];
end


figure;
histogram(m_log,20);
hold on;%xlim([0,1]);

histogram(m_perm_log,20);
legend('fitting output','fitting output permuted');

%% top/bottom 25 percentile of spillover for random smooth partition and output of model.

% dat1 = load('~/Google Drive/Presentations/Presentations/ASM_figures/modelCell_metric_fittedSU8.mat');
% dat2 = load('~/Google Drive/Presentations/Presentations/ASM_figures/modelCell_metric_smoothSU_control_fittedSU8.mat');
dat1 = load('~/Google Drive/Presentations/Presentations/ASM_figures/modelCell_metric.mat');
dat2 = load('~/Google Drive/Presentations/Presentations/ASM_figures/modelCell_metric_smoothSU_control_gaussian2.mat');

m_log=[];m_log_control=[];
for imc =1:length(dat1.mc_data)
    cell_metric_sus = dat1.mc_data(imc).metric_sus(:,2);
    cell_metric_sus1 = dat1.mc_data(imc).metric_sus(:,1);
    thr = prctile(cell_metric_sus,0);
    m_log = [m_log;mean(cell_metric_sus(cell_metric_sus>thr))];
    
    cell_metric_sus_control = dat2.mc_data(imc).metric_sus(:,2);
    cell_metric_sus_control1 = dat2.mc_data(imc).metric_sus(:,1);
    thr = prctile(cell_metric_sus_control,0);
    m_log_control = [m_log_control;mean(cell_metric_sus_control(cell_metric_sus_control>thr))];
end

figure;
histogram(m_log,20);
hold on;
histogram(m_log_control,20);
legend('simulation','control');

% take care that model, u_spatial_log and u_spatial_log_fitted are
% consistent!
plotSU_filters(u_spatial_log_gauss,model);
[metric_gauss,metric_sus_gauss] = cutting_metric(model,u_spatial_log_gauss);
cell_metric_sus = metric_sus_gauss(:,2);
cell_metric_sus1 = metric_sus_gauss(:,1);
th_log=[];th_log2=[];
for ithr =[0,25,50,75] 
thr = prctile(cell_metric_sus(:,1),ithr);
th_log = [th_log;mean(cell_metric_sus(cell_metric_sus>thr))];
th_log2 = [th_log2;mean(cell_metric_sus1(cell_metric_sus>thr))];
end
suptitle(sprintf('0: %0.02f,%0.02f ;25: %0.02f,%0.02f, 50: %0.02f,%0.02f, 75: %0.02f,%0.02f ',th_log2(1),th_log(1),th_log2(2),th_log(2),th_log2(3),th_log(3),th_log2(4),th_log(4)));

plotSU_filters(u_spatial_log_fitted,model);
[metric_fitted,metric_sus_fitted] = cutting_metric(model,u_spatial_log_fitted);
cell_metric_sus = metric_sus_fitted(:,2);
cell_metric_sus1 = metric_sus_fitted(:,1);
th_log=[];th_log2=[];
for ithr =[0,25,50,75] 
thr = prctile(cell_metric_sus(:,1),ithr);
th_log = [th_log;mean(cell_metric_sus(cell_metric_sus>thr))];
th_log2 = [th_log2;mean(cell_metric_sus1(cell_metric_sus>thr))];
end
suptitle(sprintf('0: %0.02f,%0.02f ;25: %0.02f,%0.02f, 50: %0.02f,%0.02f, 75: %0.02f,%0.02f ',th_log2(1),th_log(1),th_log2(2),th_log(2),th_log2(3),th_log(3),th_log2(4),th_log(4)));
