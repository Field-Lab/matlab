%% GLM test suite script.
% This is main script. I am trying to make it very modular, so that same
% functions can be called in other places.

addpath(genpath('../'));

%startup_tennessee
startup_bertha
%startup_bhaishahster

%% Load GLM dataset.
%load('~/Box Files Backup (not synced)/Chichilnisky Lab/ONPar_5866.mat');
load('/Volumes/Analysis/nishal/GLM_cells/nishal_glmfits/30min/7742.mat');

%% Get k, and temporal filters 
for itime=1:30
    itime
imagesc(sum(fittedGLM.cellinfo.WN_STA(:,:,:,itime),3))
colormap gray
colorbar
caxis([min(fittedGLM.cellinfo.WN_STA(:)),max(fittedGLM.cellinfo.WN_STA(:))]);
pause 
end

cell_params{1}.STAlen=30;

k=fittedGLM.linearfilters.Stimulus.Filter;
xcoords=fittedGLM.linearfilters.Stimulus.y_coord;
ycoords=fittedGLM.linearfilters.Stimulus.x_coord;


Filtdim1=size(k,1)+3;
Filtdim2=size(k,2)+3;
Filtlen = size(k,3);

stas=cell(1,1);
stas{1}=zeros(Filtdim1,Filtdim2,3,Filtlen);

% TODO : Figure out how to make 1-D STA in time into 3-D!
for itime=1:Filtlen
stas{1}(1:size(k,1),1:size(k,2),1,itime)=k(:,:,itime)'; % TODO doubt - not k/3 ??
stas{1}(1:size(k,1),1:size(k,2),2,itime)=k(:,:,itime)';
stas{1}(1:size(k,1),1:size(k,2),3,itime)=k(:,:,itime)';

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
x=sum(fittedGLM.cellinfo.WN_STA,3);
imagesc(squeeze(x(:,:,1,27)));
colormap gray
axis image
colorbar
caxis([min(x(:)),max(x(:))]);
title('Actual Spatial WN STA');

subplot(2,2,3)
plot(squeeze(response.analyse.STA(3,7,:)));
title('Simulated Temporal WN STA');

subplot(2,2,4)
plot(squeeze(sum(x(11,3,:,:),3)));
title('Actual Temporal WN STA');

%%


%% Use same code as used in experiment 
% orig_raster_sd_log=[];
% null_raster_sd_log=[];
% for imovie_run=1:5
% 
% datafile = 'load_from_cell_params';
% stas_big{1}=zeros(32,64,3,30);
% stas_big{1}(xcoords,ycoords,1,end:-1:1)=k;
% stas_big{1}(xcoords,ycoords,2,end:-1:1)=k;
% stas_big{1}(xcoords,ycoords,3,end:-1:1)=k;
% stas_big2{1}=stas_big{1}+rand(size(stas_big{1}))*0.03;
% 
% figure;
% for itime=1:30
%     itime
% imagesc(sum(stas_big2{1}(:,:,:,itime),3))
% colormap gray
% colorbar
% caxis([min(stas_big2{1}(:)),max(stas_big2{1}(:))]);
% pause(1/120); 
% end
% 
% cell_params2=struct();
% cell_params2.type_name_inp='cell';
% % cell_params.cell_list=[3888,2825,1820,4129, 5346,5671,5161,1278, 3828,3574,4036,3572, 503,560,797,1009,487,181,901]; % if type_name_inp = 'userCellList' 
% cell_params2.use_fits=2;
% cell_params2.STAlen=14;
% cell_params2.stas=stas_big2;
% 
% mov_params2=struct();
% mov_params2.mov_type='bw';
% mov_params2.movie_time=120*10;
% mov_params2.mean=0.5*255;
% mov_params2.deviation=0.48*255;
% mov_params2.scaling_loss=0.01; % a number in [0,1], fraction of values that is changed by scaling.
% 
% solver=3;
% [mov_orig,mov_modify_new]=null_space_movie2(datafile,cell_params2,mov_params2,solver);
% mov_orig=(mov_orig-127.5)/255;
% mov_modify_new =(mov_modify_new - 127.5)/255;
% 
% 
% %% Generate response
% cell_params3=cell_params;
% cell_params3{1}.stas=zeros(32,64,3,30);
% 
% for itime=1:30
%     for icol=1:3
% cell_params3{1}.stas(:,:,icol,itime)=stas_big{1}(:,:,icol,30-itime+1);
%     end
% end
% 
% 
% mov_params_orig.nTrials=50;
% mov_params_orig.mov = mov3Dto4D(mov_orig);
% mov_params_orig.movie_len = size(mov_orig,3);
% mov_params_orig.refresh = mov_params.refresh;
% response_orig=generate_response_ts(mov_params_orig,cell_params3);
% 
% mov_params_null.nTrials=50;
% mov_params_null.mov =  mov3Dto4D(mov_modify_new);
% mov_params_null.movie_len = size(mov_modify_new,3);
% mov_params_null.refresh = mov_params.refresh;
% response_null=generate_response_ts(mov_params_null,cell_params3);
% 
% 
% figure;    
% subplot(2,1,1);
% plotSpikeRaster(logical(response_orig.spksGen),'PlotType','vertline');
% title('Raster Original');
% subplot(2,1,2);
% plotSpikeRaster(logical(response_null.spksGen),'PlotType','vertline');
% title('Raster Null');
% 
% 
% %% Raster metric calculation
% metric_params.type='psth-sd';
% response_orig = calculate_raster_metric(response_orig,metric_params);
% response_null = calculate_raster_metric(response_null,metric_params);
% 
% response_log_orig{imovie_run}=response_orig;
% response_log_null{imovie_run}=response_null;
% 
% orig_raster_sd_log(imovie_run)=response_orig.psth_sd;
% null_raster_sd_log(imovie_run)=response_null.psth_sd;
% 
% end


%% Analyse actual experimental data now!
%% Actual experimental data ? 

rawMovFrames=5760;
[stim,height,width,header_size] = get_raw_movie('/Volumes/Data/2014-11-05-2/visual/18.rawMovie',rawMovFrames,1);
subtract_movies{3}=mean(stim,1);
subtract_movies{3}=mean(stim,1)*0+127.5;
movie=stim-repmat(subtract_movies{3},[rawMovFrames,1,1]);
movie=movie/255;

movStart=12*120*[0:4];
condMovies=cell(4,1);
for icond=1:4
    condMovies{icond}=movie(movStart(icond)+1:movStart(icond+1),:,:);
end

figure
for icond=1:4
    subplot(2,2,icond);
    qq=condMovies{icond}(120:12*120-120,:,:);
    hist(qq(:),20)
    title(sprintf('Condition %d',icond));
end


figure('Color','w');
for icondi=2:3
    subplot(2,1,icondi-1);
    aaa=condMovies{icondi}(240:1440-240,:,:)-condMovies{1}(240:1440-240,:,:);
    hist(aaa(:),50);
    
end
%% Condition strings
nConditions=4;
condDuration=12;
cond_str=cell(4,1);
cond_str{1}='Original';
cond_str{2}='Null for On Parasol';
cond_str{3}='Null for Off Parasol';
cond_str{4}='Null for 4 cells'
interestingConditions=[1,2,3,4];
%%
WN_datafile = '2014-11-05-2/data009/data009';
Null_datafile = '/Volumes/Analysis/2014-11-05-2/data010';
%WN_mov='/Volumes/Analysis/stimuli/white-noise-xml/RGB-10-2-0.48-11111.xml';
WN_datafile_full = '/Volumes/Analysis/2014-11-05-2/data009';

datarun=load_data(WN_datafile)
datarun=load_params(datarun)

InterestingCell_vis_id = fittedGLM.cellinfo.cid;
imov=10;
ref_cell_number=1;
[spkColl,spkCondColl]=plot_raster_script(datarun,WN_datafile,WN_datafile_full,Null_datafile,InterestingCell_vis_id,imov,ref_cell_number,nConditions,condDuration,cond_str);
     
%%
stas_big{1}=zeros(32,32,3,30);
stas_big{1}(ycoords,xcoords,1,end:-1:1)=k;
stas_big{1}(ycoords,xcoords,2,end:-1:1)=k;
stas_big{1}(ycoords,xcoords,3,end:-1:1)=k;

cell_params3=cell_params;
cell_params3{1}.stas=zeros(32,32,3,30);

for itime=1:30
    for icol=1:3
cell_params3{1}.stas(:,:,icol,itime)=stas_big{1}(:,:,icol,30-itime+1)';
    end
end
for icond=1:nConditions;
mov_use=zeros(32,32,1440);
for itime=1:size(mov_use,3)
mov_use(:,:,itime)=squeeze(condMovies{icond}(itime,:,:));
end

mov_params2.nTrials=30;
mov_params2.mov = mov3Dto4D(mov_use);
mov_params2.movie_len = size(mov_use,3);
mov_params2.refresh = mov_params.refresh;
response_pred_GLM{icond}=generate_response_ts(mov_params2,cell_params3);
end

%%
col='rkrkrkrk';
figure;    
for icond=1:4

    [x2{icond},y2{icond}]=plotSpikeRaster(logical(response_pred_GLM{icond}.spksGen),'PlotType','vertline');
    [x1{icond},y1{icond}]=plotSpikeRaster(spkCondColl(icond).spksColl,'PlotType','vertline');
end

figure;
icnt=0;
for icond=1:4
plot(x1{icond},y1{icond}+(icnt)*30,col(icnt+1));
hold on
icnt=icnt+1;
plot(x2{icond}*(20000*12)/(10*1440),y2{icond}+(icnt)*30,col(icnt+1));
hold on
icnt=icnt+1;
end

