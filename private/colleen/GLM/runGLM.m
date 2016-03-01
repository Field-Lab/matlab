% runGLM

load('/Volumes/Lab/Users/akheitman/NSEM_Home/BlockedSpikes/2012-08-09-3/NSEM_mapPRJ/organizedspikes_ONPar_5134.mat')

%2012-08-09-3
% NSEM_eye-120-3_0-3600 schemeA
load '/Volumes/Lab/Users/akheitman/NSEM_Home/Stimuli/NSEM_eye-120-3_0-3600/fitmovie_schemeA_8pix_Identity_8pix.mat';

%2012-09-27-3
% NSEM_eye-120-3_0-3600 schemeA
% load '/Volumes/Lab/Users/akheitman/NSEM_Home/Stimuli/NSEM_eye-120-3_0-3600/fitmovie_schemeA_8pix_Identity_8pix.mat';

% 2013-08-19-6 
% load '/Volumes/Lab/Users/akheitman/NSEM_Home/Stimuli/NSEM_eye-long-v2/fitmovie_schemeA_8pix_Identity_8pix.mat';

%2013-10-10-0
% load '/Volumes/Lab/Users/akheitman/NSEM_Home/Stimuli/NSEM_FEM900FF_longrast/fitmovie_schemeB_8pix_Identity_8pix.mat';



blockedspikes = organizedspikes.block.t_sp_withinblock(2:2:end);
blockedmoviecell = NSEMmovie.fitmovie.movie_byblock;

[fitspikes, fitmovie] = blocked_prep(blockedspikes, blockedmoviecell);
[STA, center] = STA_Test(fitspikes, fitmovie, 0);
[fittedGLM] = glm_fit(fitspikes, fitmovie, center);
plotfilters(fittedGLM);