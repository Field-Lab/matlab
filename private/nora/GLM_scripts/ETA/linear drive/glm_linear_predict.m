%% NB 2015-05-04
% This function takes in a GLM fit and a movie and returns the linear drive

function lcif = glm_linear_predict(fittedGLM,testmovie)
%%
% INPUTS
% fittedGLM structure
% testmovie should be in stim size x time (NO RGB!)

%% Get some info from the fittedGLM
bpf               = fittedGLM.bins_per_frame;
params.bins       = fittedGLM.bins_per_frame *size(testmovie,3);
params.frames     = size(testmovie,3);
center_coord = fittedGLM.center_coord;
frame_shifts = fittedGLM.linearfilters.Stimulus.frame_shifts;
ROI_pixels   = length(fittedGLM.linearfilters.Stimulus.x_coord) *length(fittedGLM.linearfilters.Stimulus.y_coord); 

%% Stimulus Prep
GLMType_fortest                 = fittedGLM.GLMType;
GLMType_fortest.stimfilter_mode = 'fullrank';   % treat all filters the same
[X_frame] = prep_stimcelldependentGPXV(GLMType_fortest, fittedGLM.GLMPars, testmovie,fittedGLM.inputstats,center_coord) ;
clear GLMType_fortest
GLMType = fittedGLM.GLMType;

  
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
% Nonlinearity
% if isfield(GLMType, 'lcif_nonlinearity')
%     lcif_kx_frame0 = lcif_kx_frame;
%     if strcmp(GLMType.lcif_nonlinearity.type,'piece_linear_aboutmean')
%         par = GLMType.lcif_nonlinearity.increment_to_decrement;
%         pos_mult  = (2*par) / (par + 1) ;
%         neg_mult  =      2  / (par + 1) ;        
%         pos_ind = find(lcif_kx_frame0>0);
%         neg_ind = find(lcif_kx_frame0<=0);
%         lcif_kx_frame = lcif_kx_frame0;
%         lcif_kx_frame(pos_ind) = pos_mult * (lcif_kx_frame(pos_ind));
%         lcif_kx_frame(neg_ind) = neg_mult * (lcif_kx_frame(neg_ind));
%     elseif strcmp(GLMType.lcif_nonlinearity.type,'oddfunc_powerraise_aboutmean')
%         par = GLMType.lcif_nonlinearity.scalar_raisedpower;
%         pos_ind = find(lcif_kx_frame0>0);
%         neg_ind = find(lcif_kx_frame0<=0);
%         lcif_kx_frame = lcif_kx_frame0;
%         lcif_kx_frame(pos_ind) =  (     (lcif_kx_frame(pos_ind))  .*par );
%         lcif_kx_frame(neg_ind) = -( (abs(lcif_kx_frame(neg_ind))) .*par );
%     end
% end

%% Get the drive!
lcif_kx0 = reshape( repmat(lcif_kx_frame, bpf, 1) , 1 , params.bins);
lcif_mu0 = MU * ones (1,params.bins);
lcif  = lcif_kx0; %+ lcif_mu0;

end





