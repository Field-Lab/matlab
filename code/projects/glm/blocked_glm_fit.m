tic

disp('Loading spikes and movie')
load('/Volumes/Lab/Users/akheitman/NSEM_Home/BlockedSpikes/2012-08-09-3/WN_mapPRJ/organizedspikes_ONPar_841.mat')
spikes = organizedspikes.block.t_sp_withinblock(2:2:end);  % the odd blocks are the rasters
load('/Volumes/Lab/Users/akheitman/NSEM_Home/Stimuli/BW-8-1-0.48-11111_RNG_16807/fitmovie_8pix_Identity_8pix.mat')
center = [50, 21];

%%
toc
disp('Concatenating the blocks')
[fitspikes, fitmovie] = blocked_prep(spikes, BWmovie.fitmovie.movie_byblock);

%% Just to check!
toc
STA = STA_Test(fitspikes, fitmovie, center);

%%
toc
disp('Fitting time!')
fittedGLM = glm_fit(fitspikes, fitmovie, center, 'WN_STA',STA);
plotfilters(fittedGLM);