%load('/Volumes/Lab/Users/Nora/NSEM_Home/GLMOutput_Raw/rk1_MU_PS_noCP_SUexp_init_p8IDp8fit/standardparams/NSEM_mapPRJ/2012-09-27-3/ONPar_91.mat')
load('/Volumes/Lab/Users/Nora/NSEM_Home/GLMOutput_Raw/rk1_MU_PS_noCP_SUexp_init_p8IDp8fit/standardparams/NSEM_mapPRJ/2012-08-09-3/ONPar_841.mat')

[PSTH_sim_SU, PSTH_rec] = plotrasters(fittedGLM.xvalperformance, fittedGLM, 'separate', true);

%load('/Volumes/Lab/Users/Nora/NSEM_Home/GLMOutput_Raw/rk1_MU_PS_noCP_p8IDp8/standardparams/NSEM_mapPRJ/2012-09-27-3/ONPar_91.mat')
load('/Volumes/Lab/Users/Nora/NSEM_Home/GLMOutput_Raw/rk1_MU_PS_noCP_p8IDp8/standardparams/NSEM_mapPRJ/2012-08-09-3/ONPar_841.mat')
[PSTH_sim_rec] = plotrasters(fittedGLM.xvalperformance, fittedGLM, 'separate', true);
