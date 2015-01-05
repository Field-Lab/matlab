%% Load the GLM fit
load('/Volumes/Analysis/nora/NSEM/GLM_Output/old_fits/rk1_MU_PS_noCP_p8IDp8/standardparams/NSEM_mapPRJ/2012-08-09-3/OFFPar_1471.mat');

%% Load the datarun to compare to
datarun=load_data('2012-08-09-3/data005');
datarun=load_neurons(datarun);

%% Load the movie to make predictions for
[testmovie] = loadmoviematfile(fittedGLM.cellinfo.exp_nm , 'NSEM', '8pix_Identity_8pix','testmovie');

%% Evaluate
x=GLM_predict(fittedGLM, datarun, testmovie{1,1}.matrix, 61);
plotraster(x,fittedGLM,'labels',true)
