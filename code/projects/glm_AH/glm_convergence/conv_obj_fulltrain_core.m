function [objval logprob bps] = conv_obj_fulltrain_core(GLMType,P_opt, spikes, fitmovie, inputstats, cellinfo)
% AKHEITMAN
% 2015-02-06 Compute Full objective function
logprob = [];
bps = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load up GLMParams compute some universal params
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
GLMPars           = GLMParams;
if isfield(GLMType, 'specialchange') && GLMType.specialchange
    GLMPars = GLMParams(GLMType.specialchange_name);
end

frames = size(fitmovie,3);
bins   = frames * GLMPars.bins_per_frame;
t_bin  = cellinfo.computedtstim / GLMPars.bins_per_frame; % USE THIS tstim!! %

% Perhaps we should combine this! With convolving with spikes !
bin_size      = t_bin;
if GLMType.PostSpikeFilter
    basis_params  = GLMPars.spikefilters.ps;
    ps_basis      = prep_spikefilterbasisGP(basis_params,bin_size);
end
if GLMType.CouplingFilters
    basis_params  = GLMPars.spikefilters.cp;
    cp_basis      = prep_spikefilterbasisGP(basis_params,bin_size);
end
clear bin_size basis_params

% Convolve Spike Times with appropriate basis
% Think about flushing dt out to the wrapper
% Take care of all timing in glm_execute or in glmwrap.
t_bin        = t_bin;
home_sptimes = spikes.home';
home_spbins  = ceil(home_sptimes / t_bin);
home_spbins = home_spbins(find(home_spbins < bins) );
if GLMType.PostSpikeFilter
    basis         = ps_basis';
    PS_bin        = prep_convolvespikes_basis(home_spbins,basis,bins);
end
if GLMType.CouplingFilters;
    basis = cp_basis';
    display('figure out coupling here!  CP_bin');
end
if GLMType.TonicDrive
    MU_bin = ones(1,bins);
end


%%
% PREPARE PARAMETERS
[paramind] =  prep_paramindGP(GLMType, GLMPars); 
% ORGANIZE STIMULUS COVARIATES
center_coord       = cellinfo.slave_centercoord;
WN_STA             = double(cellinfo.WN_STA);
[X_frame,X_bin]    = prep_stimcelldependentGPXV(GLMType, GLMPars, fitmovie, inputstats, center_coord, WN_STA);
clear WN_STA center_coord



%% Run through optimization .. get out pstart, fstar, eflag, output
% CONVEXT OPTIMIZATION
if GLMType.CONVEX
    glm_covariate_vec = NaN(paramind.paramcount , bins );  % make sure it crasheds if not filled out properly
    % Maybe move this inside of the stimulus preparation % 
    bpf         = GLMPars.bins_per_frame;
    shifts      = 0:bpf:(GLMPars.stimfilter.frames-1)*bpf;
    if isfield(GLMPars.stimfilter,'frames_negative')
        shifts = -(GLMPars.stimfilter.frames_negative)*bpf:bpf:(GLMPars.stimfilter.frames-1)*bpf;
    end
    
    if strcmp(GLMType.stimfilter_mode, 'fixedSP_rk1_linear')
        X_bin_shift = prep_timeshift(X_bin,shifts);
    elseif strcmp(GLMType.stimfilter_mode, 'fixedSP-ConductanceBased')
        X_bin_shift_E = prep_timeshift(X_bin,shifts);
        X_bin_shift_I = prep_timeshift(X_bin,shifts); 
        X_bin_shift = [X_bin_shift_E ; X_bin_shift_I];
        
        nonlinearity.type         = 'ConductanceBased_HardRect';
        nonlinearity.linear_index = 1:((paramind.X)-1);
        nonlinearity.excitatoryfilter_index = paramind.excitatoryfilter_index;
        nonlinearity.inhibitoryfilter_index = paramind.inhibitoryfilter_index;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if isfield(paramind, 'MU')
        glm_covariate_vec( paramind.MU , : ) = MU_bin;
    end
    if isfield(paramind, 'X')
        glm_covariate_vec( paramind.X , : ) = X_bin_shift;
    end
    if isfield(paramind, 'PS')
        glm_covariate_vec( paramind.PS , : ) = PS_bin;
    end
    if isfield(paramind, 'CP')
        glm_covariate_vec( paramind.CP , : ) = CP_bin;
    end
    
    objval = NaN(length(P_opt),1);
    for i_param = 1:length(P_opt)
        fit_param = P_opt{i_param};
        f_eval = compute_onlyf_linear(fit_param, glm_covariate_vec,home_spbins,t_bin);
        objval(i_param) = f_eval;
    end
    
    if nargin > 1
        spikecount = length(home_spbins);
        adj_factor = spikecount*log(t_bin);
        logprob_conv = adj_factor - objval;
        sprate_bin  = spikecount / bins;
        sprate_secs = sprate_bin /t_bin;
        p0         = zeros(length(fit_param),1);
        p0(1)       = log(sprate_secs);
        f_eval_null = compute_onlyf_linear(p0, glm_covariate_vec,home_spbins,t_bin);
        logprob_null = adj_factor - f_eval_null;
        bps_conv = (logprob_conv - logprob_null) / spikecount;
        
        
        logprob = logprob_conv;
        bps = bps_conv;
    end
    

end

end

function f_eval = compute_onlyf_linear(linear_params, covariates,spikebins,bin_duration)
    p = linear_params;
	COV = covariates;
	dt = bin_duration;
	spt = spikebins;
     % Find Conditional Intensity and its log
	lcif = p' * COV;
	cif  = exp(lcif);
    % Evaluate the objective function (monotonic in log-likelihood)
	f_eval = sum( lcif(spt) ) - dt * sum(cif);
    
    f_eval = -f_eval;
end