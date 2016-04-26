%% NB 2015-05-04

function [lcif] = glm_gen_signal_spatial(fittedGLM,testmovie)
%%
try
    center_coord = fittedGLM.center_coord;
catch
    center_coord = fittedGLM.cellinfo.slave_centercoord;
end
ROI_pixels   = length(fittedGLM.linearfilters.Stimulus.x_coord) *length(fittedGLM.linearfilters.Stimulus.y_coord); 

%%
GLMType_fortest                 = fittedGLM.GLMType;
GLMType_fortest.stimfilter_mode = 'fullrank';   % treat all filters the same
inputstats.mu_avgIperpix = mean(testmovie(:));
inputstats.range = max(testmovie(:))-min(testmovie(:));
[X_frame] = prep_stimcelldependentGPXV(GLMType_fortest, fittedGLM.GLMPars, testmovie,inputstats,center_coord) ;
clear GLMType_fortest
  
    %% Set up CIF Components
K = fittedGLM.linearfilters.Stimulus.space_rk1;
K  = reshape(K, [1 ROI_pixels]);
lcif = K*X_frame;

%display('binning the lcif components')
%lcif_kx0 = reshape( repmat(lcif_kx_frame, bpf, 1) , 1 , params.bins);
%lcif_mu0 = MU * ones (1,params.frames);
%lcif = lcif_kx_frame+lcif_mu0;

end





