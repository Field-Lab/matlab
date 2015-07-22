
%% Make model 
model = model_LNLN();

%% Generate test response
sz = model.gridSzX/ (8*3);


Tlen = 120*10;
movie = (randn(sz,sz,Tlen)>0)-0.5;
dt=1/120;
nTrials=30;
response = generateResp_LNLN(model,movie,dt,nTrials);
figure;
plotSpikeRaster(response~=0,'PlotType','vertline');


Tlen = 120*60*30;
movie = (randn(sz,sz,Tlen)>0)-0.5;
dt=1/120;
nTrials=1;
response = generateResp_LNLN(model,movie,dt,nTrials);
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

figure;
imagesc( (model.totalConeMap3D==0).*ssf*15 + model.totalConeMap3D)