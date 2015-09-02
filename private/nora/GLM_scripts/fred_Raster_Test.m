% cell = 'OFFPar_5161.mat';
cell = 'ONPar_3152.mat';
filename = ['/Volumes/Lab/Users/Nora/NSEM_Home/GLMOutput_Raw/rk2_MU_PS_CP_p8IDp8/standardparams/NSEM_mapPRJ/2012-08-09-3/CP_PCA/' cell];

% cell = 'OFFPar_31.mat';
% cell = 'ONPar_6858.mat';
% filename = ['/Volumes/Lab/Users/Nora/NSEM_Home/GLMOutput_Raw/rk2_MU_PS_CP_p8IDp8/standardparams/NSEM_mapPRJ/2012-09-27-3/CP_PCA/' cell];

load(filename)
plotrasters(fittedGLM.xvalperformance, fittedGLM, 'PSTH', true, 'PSTH_window_size', 200);