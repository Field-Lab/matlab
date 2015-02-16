%% For dataset 2014-11-05-2-data010. If the linear filter was estimated correctly. 

% Maybe, the response to null stimulus could be explained by the fact that
% the GLM is correct encoding model, but somehow (finite amount of data to fit it, or fundamental biases in our estimation procedure)
% make it hard for us to estimate the linear filter correctly. In this
% case, the GLM should be the correct predictor, only if the linearity was
% estimated correctly. 
% 
% Work done for cells selected for presentation by EJ.
% 
% Plan: 
% 1 .Take a GLM fitted real cell. (Keep the linear front end rank 1, rank 2 or full STA).
% 2. Generate spikes using that for 15 minutes, and 30 minutes.
% 3. Run the null space algorithm. (careful, to use the null movie (post-processing fcn, clip STA) algorithm which was used on experiment day).
% 4. Generate the rasters to original and null movie. See how much structure is there in the movie.
% 
% script taken from example_script_extended.m

%% add path to nora's folder for GLM code
location_of_git_repo='/home/vision/Nishal/matlab/';
addpath(genpath([location_of_git_repo '/private/nora']));

% Java library
javaaddpath('/Volumes/Lab/Development/vision7/Vision.app/Contents/Resources/Java/Vision.jar');
addpath(genpath('~/Nishal/matlab/private/nishal/create_act_2/'));
addpath(genpath('~/Nishal/matlab/private/nishal/Test suite/GLM/'));
addpath(genpath('~/Nishal/matlab/private/nishal/Test suite/GLM2/'));
addpath(genpath('~/Nishal/matlab/code'));

%% If want to do run code for experiment day
rmpath(genpath('~/Nishal/matlab/private/nishal/create_act_2/'));
addpath(genpath('~/Nishal/matlab/private/nishal/create_act_2_exp_repo/create_act_2_pc2014_11_05_2/'));
code_version='old';

%% Get a cell and fit GLM
cellID=1382;
% goto /matlab/code/projects/GLM/GLM_pars.m to change the 
  %fittedGLM=glm_fit_from_WNrun({cellID}, '2014-11-05-2/data009', 'RGB-10-2-0.48-11111-32x32', 1800);
% Or load a GLM fit that already exists

 load(sprintf('/Volumes/Analysis/nora/nishal_glmfits/30min/%d.mat',cellID));
 
 %%
 
sim_rast_tv1Actual=[];
sim_rast_tv2Actual=[];
sim_rast_tv1STA=[];
sim_rast_tv2STA=[];
WNstas=cell(30,1);

for imodel_run=1:30
imodel_run
    % Generate response to WN
    WNtime=120*60*30;
    WNmovie =double(rand(32,32,WNtime)>0.5)-0.5;
    x=GLM_predict(fittedGLM, WNmovie, 1);
    
    % Calculate STA 
    % WNmovie made 4 dimensional
    WNmov4D = zeros(32,32,3,WNtime);
    for iframe=1:size(WNmovie,3);
    WNmov4D(:,:,1,iframe)=WNmovie(:,:,iframe);
    WNmov4D(:,:,2,iframe)=WNmovie(:,:,iframe);
    WNmov4D(:,:,3,iframe)=WNmovie(:,:,iframe);
    end
    mov_params.mov=WNmov4D;
    
    sta_params.Filtlen=30;
    sta_params.useTrial=1;
    
    cell_params.binsPerFrame=10;
    
    response.spksGen=x.rasters.glm_sim;
    aa=repmat([1:WNtime],[10,1]);
    response.mov_frame_number=aa(:);
    
    response = calculate_sta_ts(mov_params,response,sta_params,cell_params);
    WNSTA = response.analyse.STA;
    WNstas{imodel_run}=WNSTA;
    
         figure
         for itime=1:sta_params.Filtlen
         imagesc(squeeze((WNSTA(:,:,itime)))');colormap gray
         caxis([min(WNSTA(:)),max(WNSTA(:))]);
         colorbar
         pause(1/120)
         end

    % Generate null movie from STA calculated above ? 
      mov_type_null_touse='bw';
      null_filter='STA';
      mdf_file ='/Volumes/Analysis/stimuli/white-noise-xml/BW-10-1-0.48-11111-32x32.xml';
      null_movie_compute_ts
    %  null_compute_usingactual_GLMfilter
      testmovie_filename='~/Nishal/TS_data/18.rawMovie';
      testmovie=get_rawmovie(testmovie_filename,2880);
      testmovieSTA=permute(testmovie,[2 3 1]);
       
 
% Generate null movie using actual GLM filter
      mov_type_null_touse='bw';
      null_filter='actual';
     mdf_file ='/Volumes/Analysis/stimuli/white-noise-xml/BW-10-1-0.48-11111-32x32.xml';
     null_movie_compute_ts 
      testmovie_filename='~/Nishal/TS_data/18.rawMovie';
      testmovie=get_rawmovie(testmovie_filename,2880);
      testmovieActual=permute(testmovie,[2 3 1]);
       
    % Generate rasters
% STA null
x=GLM_predict(fittedGLM, testmovieSTA, 100);
%x.rasters.recorded = x.rasters.glm_sim;
plotraster(x,fittedGLM,'labels',true,'raster_length',24)

% Actual null
y=GLM_predict(fittedGLM, testmovieActual, 100);
plotraster(y,fittedGLM,'labels',true,'raster_length',24)

% Together

x.rasters.recorded = y.rasters.glm_sim;
plotraster(x,fittedGLM,'labels',true,'raster_length',24)


 rec_rast=x.rasters.glm_sim;
calculate_psth

cond_times1 = time>=12*0+2 & time<12*1-2;
 cond_times2 = time>=12*1+2 & time<12*2-2;

sim_rast_tv1STA(imodel_run)=sqrt(var(PSTH_rec(cond_times1)));
sim_rast_tv2STA(imodel_run)=sqrt(var(PSTH_rec(cond_times2)));

rec_rast=y.rasters.glm_sim;
calculate_psth

cond_times1 = time>=12*0+2 & time<12*1-2;
 cond_times2 = time>=12*1+2 & time<12*2-2;

sim_rast_tv1Actual(imodel_run)=sqrt(var(PSTH_rec(cond_times1)));
sim_rast_tv2Actual(imodel_run)=sqrt(var(PSTH_rec(cond_times2)));

save('~/Nishal/Est_err_bw_old_code.mat','sim_rast_tv1Actual','sim_rast_tv2Actual','sim_rast_tv1STA','sim_rast_tv2STA','WNstas','fittedGLM');


end

%% figure;
load('~/Nishal/Est_err.mat');
figure('Color','w');

[X1,N1]=hist(sim_rast_tv1STA);
[X2,N2]=hist(sim_rast_tv1Actual);
[X3,N3]=hist(sim_rast_tv2STA);
[X4,N4]=hist(sim_rast_tv2Actual);

bar(N1,X1,'r');
hold on;
bar(N2,X2,'b');
hold on;
bar(N3,X3,'k');
hold on;
bar(N4,X4,'m');

legend('Original STA','Original Actual','Null STA','Null Actual')
