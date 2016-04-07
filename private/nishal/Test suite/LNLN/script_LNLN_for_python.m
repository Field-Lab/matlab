
%% Make model 
model = model_LNLN_fewSU();

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
%print(h,'-dpdf',sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/modelCellLNLN/model1/STA30.pdf'));

h=figure;
imagesc( repelem(sum(model.totalConeMap3D,3)==0,1,1,3).*ssf*20 + model.totalConeMap3D);
axis image
title('STA 20 scale');
set(gca,'xTick',[]);
set(gca,'yTick',[]);
%print(h,'-dpdf',sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/modelCellLNLN/model1/STA20.pdf'));

%% save data

 mask2 = (zeros(size(movie,1),size(movie,2)));
 mask2(1:5,1:5)=1;
 mask2=logical(mask2);
 
 maskedMov= filterMov(movie,mask2,squeeze(model.ttf));
 maskedMov2=[maskedMov;ones(1,size(maskedMov,2))];

 binnedResponses = response';

save('/Volumes/Lab/Users/bhaishahster/GITs/python_env/modelLNLN_small2_skew.mat','model','binnedResponses','maskedMov','mask2','-v7.3')

% Run python code, get fitted sub-units!
