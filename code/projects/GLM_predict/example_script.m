%% add path to nora's folder for GLM code
location_of_git_repo='/Users/Nora/Documents/MATLAB/matlab';
addpath(genpath([location_of_git_repo '/private/nora']));

%% Fit the GLM from a classification run or other white noise run
%fittedGLM=glm_fit_from_WNrun({316,7726}, '2014-11-05-2/data009_nps', 'RGB-10-2-0.48-11111-32x32', 900, '~/Desktop');
% Or load a GLM fit that already exists
 load('/Volumes/Analysis/nora/nishal_glmfits/15min/316.mat');

%% Load the rawMovie to make predictions for
testmovie_filename='/Volumes/Data/2014-11-05-2/visual/18.rawMovie';
testmovie=get_rawmovie(testmovie_filename,5760);
testmovie=permute(testmovie,[2 3 1]);

%% Load the datarun to compare (where the above movie was actually run)
datarun=load_data('2014-11-05-2/data010_from_data009_nps');
datarun=load_neurons(datarun);

%% Evaluate the GLM fit and plot the rasters
x=GLM_predict(fittedGLM, datarun, testmovie, 30);
plotraster(x,fittedGLM,'labels',true,'raster_length',12)