%% NB 2015-05-04

function [lcif] = glm_gen_signal_TEST(fittedGLM,testmovie,incl_mu)
%%

params.bins       = size(testmovie,3);
params.frames     = size(testmovie,3);
try
    center_coord = fittedGLM.center_coord;
catch
    center_coord = fittedGLM.cellinfo.slave_centercoord;
end
frame_shifts = fittedGLM.linearfilters.Stimulus.frame_shifts;
ROI_pixels   = length(fittedGLM.linearfilters.Stimulus.x_coord) *length(fittedGLM.linearfilters.Stimulus.y_coord); 

%%
GLMType_fortest                 = fittedGLM.GLMType;
GLMType_fortest.stimfilter_mode = 'fullrank';   % treat all filters the same
inputstats.mu_avgIperpix = mean(testmovie(:));
inputstats.range = max(testmovie(:))-min(testmovie(:));
[X_frame] = prep_stimcelldependentGPXV(GLMType_fortest, fittedGLM.GLMPars, testmovie,inputstats,center_coord) ;
clear GLMType_fortest
  
    %% Set up CIF Components
MU = fittedGLM.linearfilters.TonicDrive.Filter;
K  = fittedGLM.linearfilters.Stimulus.Filter;
K  = reshape(K, [ROI_pixels, length(frame_shifts)]);
KX = zeros(ROI_pixels, params.frames);
for i_pixel = 1:ROI_pixels
    X_frame_shift = prep_timeshift(X_frame(i_pixel,:),frame_shifts);
    tfilt = K(i_pixel,:);
    KX(i_pixel,:) = tfilt * X_frame_shift;
end
lcif_kx_frame = sum(KX,1);

%display('binning the lcif components')
%lcif_kx0 = reshape( repmat(lcif_kx_frame, bpf, 1) , 1 , params.bins);
if incl_mu
    lcif_mu0 = MU * ones (1,params.frames);
    lcif = lcif_kx_frame+lcif_mu0;
else
    lcif = lcif_kx_frame;
end

end





