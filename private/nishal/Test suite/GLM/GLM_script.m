%% GLM test suite script.
% This is main script. I am trying to make it very modular, so that same
% functions can be called in other places.

addpath(genpath('../'));

%startup_tennessee
startup_bertha
%startup_bhaishahster

%% Load GLM dataset.
%load('~/Box Files Backup (not synced)/Chichilnisky Lab/ONPar_5866.mat');
load('/Volumes/Analysis/nishal/GLM_cells/nishal_glmfits/30min/1382.mat');

%% Get k, and temporal filters 
WN_sta=cell(1,1);
WN_STA = double(fittedGLM.cellinfo.WN_STA - 0.5);
cell_params{1}.STAlen=30;
for itime=1:size(WN_STA,3)
WN_sta{1}(:,:,1,itime)=WN_STA (:,:,itime);
WN_sta{1}(:,:,2,itime)=WN_STA (:,:,itime);
WN_sta{1}(:,:,3,itime)=WN_STA (:,:,itime);
end
[WN_sta,mask]=clipSTAs(WN_sta,cell_params{1});

sig_stix = sum(sum(WN_sta{1},4),3)~=0;
sig_stix=sig_stix';


k=fittedGLM.linearfilters.Stimulus.Filter;
xcoords=fittedGLM.linearfilters.Stimulus.y_coord;
ycoords=fittedGLM.linearfilters.Stimulus.x_coord;
sig_stix=sig_stix(xcoords,ycoords);

mask=logical(zeros(16,16));
mask(1:13,1:13)=sig_stix;

Filtdim1=size(k,1)+3;
Filtdim2=size(k,2)+3;
Filtlen = size(k,3);

stas=cell(1,1);
stas{1}=zeros(Filtdim1,Filtdim2,3,Filtlen);

% TODO : Figure out how to make 1-D STA in time into 3-D!
for itime=1:Filtlen
stas{1}(1:13,1:13,1,itime)=k(:,:,itime)'; % TODO doubt - not k/3 ??
stas{1}(1:13,1:13,2,itime)=k(:,:,itime)';
stas{1}(1:13,1:13,3,itime)=k(:,:,itime)';

stas{1}(1:13,1:13,1,itime)=stas{1}(1:13,1:13,1,itime).*sig_stix';
stas{1}(1:13,1:13,2,itime)=stas{1}(1:13,1:13,2,itime).*sig_stix';
stas{1}(1:13,1:13,3,itime)=stas{1}(1:13,1:13,3,itime).*sig_stix';
end




postSpikeFilter = fittedGLM.linearfilters.PostSpike.Filter;
tonicDrive = fittedGLM.linearfilters.TonicDrive.Filter;
stas{1}(:,:,:,16:end)=0;
cell_params=cell(1,1);
cell_params{1}.stas=stas{1};
cell_params{1}.postSpikeFilter=postSpikeFilter;
cell_params{1}.tonicDrive=tonicDrive;
cell_params{1}.binsPerFrame=10;


% figures
figure('Color','w');
for itime=1:Filtlen
imagesc(sum(stas{1}(:,:,:,itime),3));
caxis([min(stas{1}(:)),max(stas{1}(:))])
colorbar
colormap gray
pause(1/120)
end

% Plot filters
figure('Color','w');
subplot(4,1,1);
plot(cell_params{1}.postSpikeFilter);
title('Post spike filter');

subplot(4,1,2);
plot(squeeze(cell_params{1}.stas(8,8,1,:)));
title('STA temporal filter');

subplot(4,1,[3,4]);
imagesc(sum(cell_params{1}.stas(:,:,:,4),3));
axis image
colormap gray
colorbar
title('STA spatial filter');


%% Generate small stimulus to generate rasters
mov_params.type='bw';
mov_params.movie_spec = '/Volumes/Analysis/stimuli/white-noise-xml/BW-8-1-0.48-11111.xml';
mov_params.movie_len =15; % in seconds
mov_params.refresh=1000/120;
mov_params = generate_movie_ts(mov_params);

mov_params.mov=mov_params.mov(1:Filtdim1,1:Filtdim2,:,:);
figure;
for itime=1:10
imagesc(mov_params.mov(:,:,itime));
colormap gray
pause(1/120)
end

% Generate response to stimulus - use k, temporal filters and movie

mov_params.nTrials=50;
response=generate_response_ts(mov_params,cell_params);
figure('Color','w');    
plotSpikeRaster(logical(response.spksGen),'PlotType','vertline');
title('Raster');

%% Generate long stimulus to calculate STA.
mov_params.type='bw';
mov_params.movie_spec = '/Volumes/Analysis/stimuli/white-noise-xml/BW-8-1-0.48-11111.xml';
mov_params.movie_len =60*30; % in seconds
mov_params.refresh=1000/120;
mov_params = generate_movie_ts(mov_params);

mov_params.mov=mov_params.mov(1:Filtdim1,1:Filtdim2,:,:);
% figure;
% for itime=1:10
% imagesc(mov_params.mov(:,:,itime));
% colormap gray
% pause(1/120)
% end

% Generate response
mov_params.nTrials=1;
response=generate_response_ts(mov_params,cell_params);

% Calculate STA
sta_params.Filtlen=40;
sta_params.useTrial=1;
response = calculate_sta_ts(mov_params,response,sta_params,cell_params{1})
%%
   
figure('Color','w');

subplot(2,2,1);
imagesc(response.analyse.STA(:,:,23));
colormap gray
axis image
colorbar
caxis([min(response.analyse.STA(:)),max(response.analyse.STA(:))]);
title('Simulated Spatial WN STA');

subplot(2,2,2)
x=WN_STA(ycoords,xcoords,:);
imagesc(x(end:-1:1,end:-1:1,5));
colormap gray
axis image
colorbar
caxis([min(x(:)),max(x(:))]);
title('Actual Spatial WN STA');

subplot(2,2,3)
plot(squeeze(response.analyse.STA(8,8,:)));
title('Simulated Temporal WN STA');

subplot(2,2,4)
plot(squeeze(x(8,8,:)));
title('Actual Temporal WN STA');

%% Generate null stimulus - use STA and movie ? 

null_mov_params.movie_idx=2;
null_mov_params.movie_len=15; % in seconds;
null_mov_params.deviation=0.36*255;
null_mov_params.scaling_loss=0.03;
null_mov_params.mov_type='bw';
null_mov_params.mean=0.5*255;
null_mov_params.totalMaskAccept=mask;

null_cell_params.STA=(mean(cell_params{1}.stas,3));
null_cell_params.Filtdim1=size(null_cell_params.STA,1);
null_cell_params.Filtdim2=size(null_cell_params.STA,2);
null_cell_params.Filtlen=size(null_cell_params.STA,4);
null_cell_params.sta_type=5;

[mov_orig,mov_modify_new]= generate_null_movie_ts(null_mov_params,null_cell_params);
%% Generate response to null stimulus

mov_params_orig.nTrials=50;
mov_params_orig.mov = mov3Dto4D(mov_orig);
mov_params_orig.movie_len = size(mov_orig,3);
mov_params_orig.refresh = mov_params.refresh;
response_orig=generate_response_ts(mov_params_orig,cell_params);

mov_params_null.nTrials=50;
mov_params_null.mov =  mov3Dto4D(mov_modify_new);
mov_params_null.movie_len = size(mov_modify_new,3);
mov_params_null.refresh = mov_params.refresh;
response_null=generate_response_ts(mov_params_null,cell_params);


figure;    
subplot(2,1,1);
plotSpikeRaster(logical(response_orig.spksGen),'PlotType','vertline');
title('Raster Original');
subplot(2,1,2);
plotSpikeRaster(logical(response_null.spksGen),'PlotType','vertline');
title('Raster Null');

