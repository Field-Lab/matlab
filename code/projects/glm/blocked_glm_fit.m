load('/Volumes/Lab/Users/akheitman/NSEM_Home/BlockedSpikes/2012-08-09-3/WN_mapPRJ/organizedspikes_ONPar_841.mat')
fitspikes = organizedspikes.block.t_sp_withinblock(2:2:end);  % the odd blocks are the rasters
load('/Volumes/Lab/Users/akheitman/NSEM_Home/Stimuli/BW-8-1-0.48-11111_RNG_16807/fitmovie_8pix_Identity_8pix.mat')
% 
% stop = 0;
% while stop == 0
%     i = 1;
%     try 
%         fitmovie{i} = NSEMmovie.fitmovie.movie_byblock{i}.matrix;
%         i = i+1;
%     catch
%         disp(['Movie is ' i ' blocks long'])
%         stop = 1;
%     end
% end

center = [20 49];


%%
fittedGLM = glm_fit(fitspikes, BWmovie.fitmovie.movie_byblock, center, 'WN_STA',STA);
plotfilters(fittedGLM);