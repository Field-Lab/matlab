


% NEW LOOK on 2015-01-
function [f grad Hess log_cif]= glm_nonconvex_optimizationfunction(params,stimfilter_mode,paramind,convex_covariates,X_frame,frame_shifts, bpf, spikebins,bin_duration,lcif_nonlinearity)
 %%% PURPOSE %%%
% HANDLE NON_CONVEX PARAMETER USAGE
% PRIMARILY RK1 RK2 STIMULUS FILTERS

% Compute the Objective Function being optimized (f)
% Compute the grad/hess as well for the optimizer
% Monotonically related to  negative of log_cif
% log_cif:= log of the conditional intensity function


%%% NOTES %%%%
% Row indexes the parameters and corresponding covariate
% Column indexes time

%%% INPUTS  %%%
% Params: glm parameters to be optimized: convex and non-convex
% Convex_Covariates: time dependent input, 
% SpikeBins: spike time in bins of cell
% Bin_Duration: duration of each bin in seconds 
% Stimfilter_mode: string determining how filter operates 
% lcif_nonlinearity: optional argument for the case of non-linearity in
%    post-convolution
% BPF: bins per frame
% X_frame: frame by frame stimulus
% frame_shifts: corresponds to how stimulus filter acts



%%% Older Comments %%%
%%% From 2014-12-15
%{
%AKHeitman 2014-05-21
% takes care of the non-convex term
% then plug into glm_convex_optimizationfunction
% this also needs to get commented better
% WORKS! 


% restarted: AKHEITMAN 2014-12-15
% CALLS
% glm_convex_optimizationfunction


%%% PURPOSE %%
% Unwrap the non-convexities, then plug into the convex opt function
% Must be inside the optimizer
% lcif_nonlinearity is new optional argument used for nonlinearties in how
% the lcif is computed
% Only call glm_convex_optimizationfunction
% rk1 works  2014-05-21
% Need to check on rk2 !!
% WORKS FOR WN, MAYBE NOT FOR NSEM
%}

% AKHEITMAN 2014-12-15  up to version_0
% 2015-01-16 starting to rethink to version_1

% Find basic dimensions
frames = size(X_frame,2);
pixels = size(X_frame,1);
bins   = frames*bpf;



% ASSIGN CONVEX COVARIATES % 
covariate_vec               = NaN(length(params) ,bins);
convex_ind                  = paramind.convParams_ind;
covariate_vec(convex_ind,:) = convex_covariates;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% REMAKE THE COVARIATE VECTORS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(stimfilter_mode, 'rk1') || strcmp(stimfilter_mode, 'rk2') || strcmp(stimfilter_mode, 'rk2-ConductanceBased')
    
    % PARAMS TO GET INCORPORTED IN COVARIATE VEC 
    TimeFilter  = params(paramind.time1);
    SpaceFilter = params(paramind.space1);
    
    % FIND SPATIAL FILTER COVARIATE VEC (USING TIMEFILTER)
    if min(frame_shifts) == 0;
        convolvingFilter = (TimeFilter);        
    elseif min(frame_shifts) > 0
        padzeros         = zeros(min(frame_shifts),1);
        convolvingFilter = [padzeros ; (TimeFilter)]; 
    else
        error('frame_shifts should be >=0')
    end
    timeconvStim = zeros( pixels , (frames+length(convolvingFilter)-1) );
    for i_row = 1:size(X_frame,1)
        timeconvStim(i_row,:) = conv(X_frame(i_row,:) , convolvingFilter);
    end
    timeconvStim = timeconvStim(:,1:frames);
    bins   = bpf * frames;
    A      = repmat(timeconvStim, [ bpf,1]); 
    spatial_covariatevec  = reshape(A, [pixels, bins]);
    
    % FIND TEMPORAL FILTER COVARIATE VEC (USING SPATIAL FILTER)
    spaceconvStim           = SpaceFilter' * X_frame;
    B                       = repmat(spaceconvStim, [ bpf,1]); 
    spaceconvStim_bin       = reshape(B, [1, bins]);
    
    bin_shifts              = bpf *frame_shifts;
    temporal_covariatevec   = prep_timeshift(spaceconvStim_bin,bin_shifts);
    
    
    % STIMULUS ADDED TWICE / DIVIDE EACH COVARIATE VEC BY TWO
    covariate_vec(paramind.time1 ,:) =  .5*temporal_covariatevec;
    covariate_vec(paramind.space1,:) =  .5*spatial_covariatevec;
end
%%% RK2 SECOND FILTER DO THE SAME %%%
if  strcmp(stimfilter_mode, 'rk2') || strcmp(stimfilter_mode, 'rk2-ConductanceBased')% same as rk1, just throwing in the rk2 term
        % PARAMS TO GET INCORPORTED IN COVARIATE VEC 
    TimeFilter  = params(paramind.time2);
    SpaceFilter = params(paramind.space2);
    
    % FIND SPATIAL FILTER COVARIATE VEC (USING TIMEFILTER)
    if min(frame_shifts) == 0;
        convolvingFilter = (TimeFilter);        
    elseif min(frame_shifts) > 0
        padzeros         = zeros(min(frame_shifts),1);
        convolvingFilter = [padzeros ; (TimeFilter)]; 
    else
        error('frame_shifts should be >=0')
    end
    timeconvStim = zeros( pixels , (frames+length(convolvingFilter)-1) );
    for i_row = 1:size(X_frame,1)
        timeconvStim(i_row,:) = conv(X_frame(i_row,:) , convolvingFilter);
    end
    timeconvStim = timeconvStim(:,1:frames);
    bins   = bpf * frames;
    A      = repmat(timeconvStim, [ bpf,1]); 
    spatial_covariatevec  = reshape(A, [pixels, bins]);
    
    % FIND TEMPORAL FILTER COVARIATE VEC (USING SPATIAL FILTER)
    spaceconvStim           = SpaceFilter' * X_frame;
    B                       = repmat(spaceconvStim, [ bpf,1]); 
    spaceconvStim_bin       = reshape(B, [1, bins]);
    
    bin_shifts              = bpf *frame_shifts;
    temporal_covariatevec   = prep_timeshift(spaceconvStim_bin,bin_shifts);
    
    
    % STIMULUS ADDED TWICE / DIVIDE EACH COVARIATE VEC BY TWO
    covariate_vec(paramind.time2 ,:) =  .5*temporal_covariatevec;
    covariate_vec(paramind.space2,:) =  .5*spatial_covariatevec;
end
%%

% Now that we have the covariate vector we can evaluate the function
if exist('lcif_nonlinearity','var')
    %display('lcif_nonlinearity')
    [f grad Hess_Base log_cif COV_NL]= glm_convex_optimizationfunction_withNL(params,covariate_vec,spikebins,bin_duration,lcif_nonlinearity);
else
    [f grad Hess_Base log_cif]= glm_convex_optimizationfunction(params,covariate_vec,spikebins,bin_duration);
end
%% NEED TO FIGURE OUT THE HESS CORRECTIONS SOME TIME .. works without
Hess = Hess_Base;
% Address non-convex terms.. dependent on filtert type.


%%% HESSIAN CORRECTION %%%

display('### No Hess Correction')

%{
if strcmp(stimfilter_mode, 'rk1') || strcmp(stimfilter_mode, 'rk2') || strcmp(stimfilter_mode, 'rk2-ConductanceBased')
    display('%%%%%%%NEWHessCorrection%%%%%%')
    [Hess_Correction] = glm_nonconvex_HessianCorrection(log_cif,spikebins,X_frame,frame_shifts,bpf,bin_duration);
    Hess = Hess_Base;
    
    Hess(paramind.space1,paramind.time1 ) = Hess(paramind.space1,paramind.time1 ) + Hess_Correction ;
    Hess(paramind.time1 ,paramind.space1) = Hess(paramind.time1 ,paramind.space1) + Hess_Correction';
    if strcmp(stimfilter_mode, 'rk2') || strcmp(stimfilter_mode, 'rk2-ConductanceBased')
        Hess(paramind.space2,paramind.time2 ) = Hess(paramind.space2,paramind.time2 ) + Hess_Correction ;
        Hess(paramind.time2 ,paramind.space2) = Hess(paramind.time2 ,paramind.space2) + Hess_Correction';
    end 
end
%}
    
    
    
    
    
end

