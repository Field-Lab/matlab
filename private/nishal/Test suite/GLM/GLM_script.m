%% GLM test suite script.
% This is main script. I am trying to make it very modular, so that same
% functions can be called in other places.

addpath(genpath('../'));

startup_tennessee

%startup_bertha
%startup_bhaishahster

%% Load GLM dataset.
%load('~/Box Files Backup (not synced)/Chichilnisky Lab/ONPar_5866.mat');
load('/Volumes/Analysis/nishal/GLM_cells/ONPar_5866.mat');

%% Get k, and temporal filters 
k=fittedGLM.linearfilters.Stimulus.Filter;
xcoords=fittedGLM.linearfilters.Stimulus.y_coord;
ycoords=fittedGLM.linearfilters.Stimulus.x_coord;


Filtdim1=size(fittedGLM.cellinfo.WN_STA,2);
Filtdim2=size(fittedGLM.cellinfo.WN_STA,1);
Filtlen = size(fittedGLM.cellinfo.WN_STA,3);

stas=cell(1,1);
stas{1}=zeros(Filtdim1,Filtdim2,3,Filtlen);

% TODO : Figure out how to make 1-D STA in time into 3-D!
for itime=1:Filtlen
stas{1}(xcoords,ycoords,1,itime)=k(:,:,itime)'/3;
stas{1}(xcoords,ycoords,2,itime)=k(:,:,itime)'/3;
stas{1}(xcoords,ycoords,3,itime)=k(:,:,itime)'/3;
end
cell_params{1}.STAlen=30;
stas=clipSTAs(stas,cell_params{1});


postSpikeFilter = fittedGLM.linearfilters.PostSpike.Filter;
tonicDrive = fittedGLM.linearfilters.TonicDrive.Filter;

cell_params=cell(1,1);
cell_params{1}.stas=stas{1};
cell_params{1}.postSpikeFilter=postSpikeFilter;
cell_params{1}.tonicDrive=tonicDrive;
cell_params{1}.binsPerFrame=10;

% figures
figure;
for itime=1:Filtlen
imagesc(sum(stas{1}(:,:,:,itime),3));
caxis([min(stas{1}(:)),max(stas{1}(:))])
colorbar
colormap gray
pause
end

%% Get Stimulus
mov_params.type='bw';
mov_params.movie_spec = '/Volumes/Analysis/stimuli/white-noise-xml/RGB-8-1-0.48-11111.xml';
mov_params.movie_len = 30*15; % in seconds
mov_params.refresh=1000/120;
mov_params = generate_movie_ts(mov_params);

figure;
for itime=1:10
imagesc(mov_params.mov(:,:,itime));
pause
end

%% Generate response to stimulus - use k, temporal filters and movie
response=generate_response_ts(mov_params,cell_params);

%% Calculate STA, use response and stimulus
sta_params.Filtlen=30;
response = calculate_sta_ts(mov_params,response,sta_params)
%% Generate null stimulus - use STA and movie ? 

%% Generate response to null stimulus
