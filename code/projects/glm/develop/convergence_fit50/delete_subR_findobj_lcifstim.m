function [lcif_stim, obj_val] = subR_findobj_lcifstim(lcif_external, pstar, GLMType, GLMPars, fitspikes, fitmovie, inputstats, glm_cellinfo, t_bin,bins)
% AKHEITMAN 2015-07-15
% subRoutine which will get optimized for inding input NL

% Part 1: Find Stim driven lcif (with input
% Part 3: Add in external lcif and find obj_val 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part 1: Find Stim driven lcif (with input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[paramind] =  prep_paramindGP(GLMType, GLMPars); 
center_coord       = glm_cellinfo.slave_centercoord;
WN_STA             = double(glm_cellinfo.WN_STA);
[X_frame,X_bin]    = prep_stimcelldependentGPXV(GLMType, GLMPars, fitmovie, inputstats, center_coord, WN_STA);
clear WN_STA center_coord



if GLMType.CONVEX
   
    bpf         = GLMPars.bins_per_frame;
    shifts      = 0:bpf:(GLMPars.stimfilter.frames-1)*bpf;
    if isfield(GLMPars.stimfilter,'frames_negative')
        shifts = -(GLMPars.stimfilter.frames_negative)*bpf:bpf:(GLMPars.stimfilter.frames-1)*bpf;
    end
    X_bin_shift = prep_timeshift(X_bin,shifts);
    pstar = pstar';
    lcif_stim = pstar(paramind.X) *X_bin_shift;
end


if ~GLMType.CONVEX
    % Maybe move this inside of the stimulus preparation % 
    bpf         = GLMPars.bins_per_frame;
    pstar = pstar;
    
    
    frame_shifts = 0:1:(GLMPars.stimfilter.frames-1);
    if isfield(GLMPars.stimfilter,'frames_negative')
        frame_shifts = -(GLMPars.stimfilter.frames_negative):1:(GLMPars.stimfilter.frames-1)*1;
    end
    pixels  = size(X_frame,1);
    frames  = size(X_frame,2);
    
    
    % PARAMS TO GET INCORPORTED IN COVARIATE VEC 
    TimeFilter  = pstar(paramind.time1);
    SpaceFilter = pstar(paramind.space1);
    
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
    stim_covariate_vec =  [.5*spatial_covariatevec; .5*temporal_covariatevec];    
    lcif_stim = (pstar(paramind.X))'  * stim_covariate_vec;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part 2: Add in external lcif and find obj_val
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargout > 1
    t_bin        = t_bin;
    home_sptimes = fitspikes.home';
    home_spbins  = ceil(home_sptimes / t_bin);
    home_spbins  = home_spbins(find(home_spbins < bins) );
    total_lcif = lcif_stim + lcif_external;
    cif        = exp(total_lcif);
    obj_val    = -(sum( total_lcif(home_spbins) ) - t_bin * sum(cif));
end

end