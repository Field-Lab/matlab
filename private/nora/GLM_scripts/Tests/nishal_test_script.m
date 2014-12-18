%% Loading
load('/Volumes/Analysis/nora/nishal_glmfits/15min/316.mat');

%%
datarun=load_data('2014-11-05-2/data010_from_data009_nps');
datarun=load_neurons(datarun);

%%
testmovie_filename='/Volumes/Data/2014-11-05-2/visual/18.rawMovie';
testmovie=get_rawmovie(testmovie_filename,5760);
testmovie=permute(testmovie,[2 3 1]);

%% Eval
x=nishal_test(fittedGLM, datarun, testmovie, 30);
plotraster(x,fittedGLM,1,0,1)
