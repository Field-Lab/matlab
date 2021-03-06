%AKHeitman 2014-05-21
% takes care of the non-convex term
% then plug into glm_convex_optimizationfunction
% this also needs to get commented better
% WORKS! 

% Only call glm_convex_optimizationfunction

% rk1 works  2014-05-21
% Need to check on rk2 !!
function [f grad Hess log_cif]= glm_nonconvex_optimizationfunction(params,stimfilter_mode,paramind,convex_covariates,X_frame,frame_shifts, bpf, spikebins,bin_duration)


% Find basic dimensions
frames = size(X_frame,2);
pixels = size(X_frame,1);
bins   = frames*bpf;

% Set up the final covariate vec. 
% Plug in the convex terms that don't need to be processed
covariate_vec               = NaN(length(params) ,bins);
convex_ind                  = paramind.convParams_ind;
covariate_vec(convex_ind,:) = convex_covariates;

% Address non-convex terms.. dependent on filtert tpe.
if strcmp(stimfilter_mode, 'rk1') || strcmp(stimfilter_mode, 'rk2')
    
    % Params already predict filters which need to multiply with movie 
    TimeFilter  = params(paramind.time1);
    SpaceFilter = params(paramind.space1);
    
    % Use TimeFilter to find Spatial Covariate Vec
    if min(frame_shifts) == 0;
        convolvingFilter = (TimeFilter);        
    elseif min(frame_shifts) > 0
        padzeros         = zeros(min(frame_shifts),1);
        convolvingFilter = [padzeros ; (TimeFilter)]; 
    else
        error('frame_shifts should be >=0')
    end
    timeconvStim = zeros( pixels , (frames +length(convolvingFilter)-1) );
    for i_row = 1:size(X_frame,1)
        timeconvStim(i_row,:) = conv(X_frame(i_row,:) , convolvingFilter);
    end
    timeconvStim = timeconvStim(:,1:frames);
    bins   = bpf * frames;
    A      = repmat(timeconvStim, [ bpf,1]); 
    spatial_covariatevec  = reshape(A, [pixels, bins]);
    
    % Use SpaceFilter to find Temporal Covariate Vec
    spaceconvStim           = SpaceFilter' * X_frame;
    B                       = repmat(spaceconvStim, [ bpf,1]); 
    spaceconvStim_bin       = reshape(B, [1, bins]);
    
    bin_shifts              = bpf *frame_shifts;
    temporal_covariatevec   = prep_timeshift(spaceconvStim_bin,bin_shifts);
    
    covariate_vec(paramind.time1 ,:) = temporal_covariatevec;
    covariate_vec(paramind.space1,:) =  spatial_covariatevec;
end

if  strcmp(stimfilter_mode, 'rk2') % same as rk1, just throwing in the rk2 term
    % TimeFilter 
    TimeFilter  = params(paramind.time2);
    SpaceFilter = params(paramind.space2);
    
    
    % Use TimeFilter to find Spatial Covariate Vec
    if min(frame_shifts) == 0;
        convolvingFilter = (TimeFilter);        
    elseif min(frame_shifts) > 0
        padzeros         = zeros(min(frame_shifts),1);
        convolvingFilter = [padzeros ; (TimeFilter)]; 
    else
        error('frame_shifts should be >=0')
    end
    
    timeconvStim = zeros( pixels , (frames +length(convolvingFilter)-1) );
    for i_row = 1:size(X_frame,1)
        timeconvStim(i_row,:) = conv(X_frame(i_row,:) , convolvingFilter);
    end
    timeconvStim = timeconvStim(:,1:frames);
    bins   = bpf * frames;
    A      = repmat(timeconvStim, [ bpf,1]); 
    spatial_covariatevec  = reshape(A, [pixels, bins]);
    
    % Use SpaceFilter to find Spatial Covariate Vec
    spaceconvStim           = SpaceFilter' * X_frame;
    B                       = repmat(spaceconvStim, [ bpf,1]); 
    spaceconvStim_bin       = reshape(B, [1, bins]);
    
    bin_shifts              = bpf *frame_shifts;
    temporal_covariatevec   = prep_timeshift(spaceconvStim_bin,bin_shifts);
    
    covariate_vec(paramind.time2 ,:) = temporal_covariatevec;
    covariate_vec(paramind.space2,:) =  spatial_covariatevec;
end

% Now that we have the covariate vector we can evaluate the function
[f grad Hess log_cif]= glm_convex_optimizationfunction(params,covariate_vec,spikebins,bin_duration);
end