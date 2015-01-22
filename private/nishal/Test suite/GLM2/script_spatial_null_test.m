%% add path to nora's folder for GLM code
location_of_git_repo='/home/vision/Nishal/matlab/';
addpath(genpath([location_of_git_repo '/private/nora']));

% Java library
javaaddpath('/Volumes/Lab/Development/vision7/Vision.app/Contents/Resources/Java/Vision.jar');
addpath(genpath('../../create_act2'));
addpath(genpath('../GLM'));

%%
cellID=3152;

load(sprintf('/Volumes/Analysis/nora/nishal_glmfits/15min/%d.mat',cellID));


    %% Generate response to WN
    WNtime=120*60*30;
    WNmovie =double(rand(32,32,WNtime)>0.5)-0.5;
    x=GLM_predict(fittedGLM, WNmovie, 1);
    
    %% Calculate STA 
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
%          figure
%          for itime=1:sta_params.Filtlen
%          imagesc(squeeze((WNSTA(:,:,itime)))');colormap gray
%          caxis([min(WNSTA(:)),max(WNSTA(:))]);
%          colorbar
%          pause(1/120)
%          end

    %% Generate null movie from STA calculated above ? 
     % null_movie_compute_ts_spatial
      null_movie_compute_ts
      testmovie_filename='~/Nishal/TS_data/18.rawMovie';
      testmovie=get_rawmovie(testmovie_filename,2880);
      testmovie=permute(testmovie,[2 3 1]);
       
    %% Generate rasters
x=GLM_predict(fittedGLM, testmovie, 30);
x.rasters.recorded = x.rasters.glm_sim;
plotraster(x,fittedGLM,'labels',true,'raster_length',24)

