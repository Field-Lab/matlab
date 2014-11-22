%% GLM test suite script.
% This is main script. I am trying to make it very modular, so that same
% functions can be called in other places.

addpath(genpath('../'));
startup_null_analyse_tenessee
startup_null_bhaishahster
%% Load GLM dataset.
%load('~/Box Files Backup (not synced)/Chichilnisky Lab/ONPar_5866.mat');
load('/Volumes/Analysis/nishal/GLM_cells/ONPar_5866.mat');

%% Get k, and temporal filters 
k=fittedGLM.linearfilters.Stimulus.Filter;
xcoords=fittedGLM.linearfilters.Stimulus.x_coord;
ycoords=fittedGLM.linearfilters.Stimulus.y_coord;


Filtdim1=size(fittedGLM.cellinfo.WN_STA,1);
Filtdim2=size(fittedGLM.cellinfo.WN_STA,2);
Filtlen = size(fittedGLM.cellinfo.WN_STA,3);

stas=cell(1,1);
stas{1}=zeros(Filtdim1,Filtdim2,3,Filtlen);

% TODO : Figure out how to make 1-D STA in time into 3-D!
for itime=1:Filtlen
stas{1}(xcoords,ycoords,1,itime)=k(:,:,itime)/3;
stas{1}(xcoords,ycoords,1,itime)=k(:,:,itime)/3;
stas{1}(xcoords,ycoords,1,itime)=k(:,:,itime)/3;
end

postSpikeFilter = fittedGLM.linearfilters.PostSpike.Filter;
tonicDrive = fittedGLM.linearfilters.TonicDrive.Filter;


%% Get Stimulus
mov_params.type='bw';
mov_params.movie_spec = '/Volumes/Analysis/stimuli/white-noise-xml/RGB-8-1-0.48-11111.xml';
mov_params.movie_len = 30*60; % in seconds
mov_params.refresh=1000/120;
mov = generate_movie_ts(mov_params);


%% Generate response to stimulus - use k, temporal filters and movie

%% Calculate STA, use response and stimulus

%% Generate null stimulus - use STA and movie ? 

%% Generate response to null stimulus
