function xvalperformance = GLM_predict(fittedGLM, testmovie, trials, datarun, start_times)
% GLM_predict(fittedGLM, datarun, testmovie, trials, optional: start_times) 
% 
% fittedGLM, the GLM structure from AH's GLM code
% testmovie, as a matfile, must be same dimensions as the fit movie
% trials, the number of trials
% datarun, optional, which must have the neurons loaded if you want the
%   recored rasters
% start_times, optional, the start time of the trials in seconds. 
%   IF each trials starts with a new trigger (wait-trigger in lisp code), 
%   you can leave this out, and the code automatically detects new trials by 
%   looking for unusually sized time intervals between triggers

%% Set some parameters


% Determine if we are getting data rasters, and if so, get
% the trial start times
if nargin==3
    predict_only=true;
elseif nargin==4
    predict_only=false;
    trigger_eps=0.2;
    trial_starts=datarun.triggers([true; abs(diff(datarun.triggers)-median(diff(datarun.triggers))) > trigger_eps ]);
elseif nargin==5
    predict_only=false;
    trial_starts=start_times;
end

% If we want data rasters to compare, then get the spikes and check trial start times 
if ~predict_only
    spikes=datarun.spikes{datarun.cell_ids==fittedGLM.cellinfo.cid};
    if length(trial_starts)<trials
        error('Did not find all trial start times')
    elseif length(trial_starts)>trials
        trial_starts=trial_starts(1:2:end);
    end
    trial_starts(end+1)=Inf;
end

testframes=size(testmovie,3);

% get some parameters from fitted GLM
bpf = fittedGLM.bins_per_frame;
params.bindur = fittedGLM.t_bin;
params.bins = fittedGLM.bins_per_frame * testframes;
params.trials = trials;
params.frames = testframes;
params.testdur_seconds = params.bindur * params.bins ;
center_coord = fittedGLM.cellinfo.slave_centercoord;
frame_shifts = fittedGLM.linearfilters.Stimulus.frame_shifts;
ROI_pixels   = length(fittedGLM.linearfilters.Stimulus.x_coord) *length(fittedGLM.linearfilters.Stimulus.y_coord) ;


%% Make the data rasters

if ~predict_only
    logicalspike = zeros( params.trials , params.bins) ;
    for i_blk = 1 : params.trials
        sptimes = spikes(spikes>trial_starts(i_blk) & spikes<trial_starts(i_blk+1))-trial_starts(i_blk);
        for i_sp = 1:length(sptimes)
            spt = sptimes(i_sp);
            binnumber = ceil(spt / params.bindur );
            try
                logicalspike( i_blk, binnumber )  =  logicalspike( i_blk,binnumber ) + 1;
            catch
            end
        end
    end
    clear i_blk spt sptimes
    
end

%% Prepare the stimulus

GLMType_fortest                 = fittedGLM.GLMType;
GLMType_fortest.stimfilter_mode = 'fullrank';   % treat all filters the same

[X_frame0 ] = prep_stimcelldependentGP(GLMType_fortest, fittedGLM.GLMPars, testmovie,center_coord) ;
X_frame     = X_frame0(:,1:testframes);

clear GLMType_fortest

% NBCoupling stuffs 06-25-2014
%{
if GLMType.CouplingFilters
    for pair=1:12
	if neighborspikes{pair} == 0
        pairspike{pair} = zeros( length(params.evalblocks) , params.bins) ;
	else
        pairspike{pair} = zeros( length(params.evalblocks) , params.bins) ;
        for i_blk = 1 : length(params.evalblocks)
            blknum = params.evalblocks(i_blk);
            sptimes = neighborspikes{pair}.block.t_sp_withinblock{blknum} - SPars.fittest_skipseconds;
            sptimes = sptimes(find(sptimes > 0 ) );
            for i_sp = 1:length(sptimes)
                spt = sptimes(i_sp);
                binnumber = ceil(spt / params.bindur );
                pairspike{pair}( i_blk, binnumber )  =  pairspike{pair}( i_blk,binnumber ) + 1;
            end
        end
	end
        clear i_blk spt sptimes
    end
end
%}

%% Set up CIF Components for prediction rasters

MU = fittedGLM.linearfilters.TonicDrive.Filter;

if fittedGLM.GLMType.PostSpikeFilter
    PS = fittedGLM.linearfilters.PostSpike.Filter;
end

%{
% NBCoupling 06-23-2014
if GLMType.CouplingFilters
    CP = fittedGLM.linearfilters.Coupling.Filter;
end
% NBcoupling
%}

K  = fittedGLM.linearfilters.Stimulus.Filter;
K  = reshape(K, [ROI_pixels, length(frame_shifts)]);

% HUGE HACK AKHeitman 2014-10-21
% rk1 filters are misscaled... too hard to dig out
% rk1 filters are fit fine
% this is confirmed to be the correct factor though!
if strcmp(fittedGLM.GLMType.stimfilter_mode, 'rk1')
    K = 2*K;
end

KX = zeros(ROI_pixels, params.frames);
for i_pixel = 1:ROI_pixels
    X_frame_shift = prep_timeshift(X_frame(i_pixel,:),frame_shifts);
    tfilt = K(i_pixel,:);
    KX(i_pixel,:) = tfilt * X_frame_shift;
end
lcif_kx_frame = sum(KX,1);

%%
%display('binning the lcif components')
lcif_kx0 = reshape( repmat(lcif_kx_frame, bpf, 1) , 1 , params.bins);
lcif_mu0 = MU * ones (1,params.bins);
lcif_mu = repmat(lcif_mu0 , params.trials, 1);
lcif_kx = repmat(lcif_kx0 , params.trials, 1);
clear sbpf;
lcif = lcif_mu + lcif_kx;
if fittedGLM.GLMType.PostSpikeFilter && ~predict_only
    lcif_ps = fastconvAH(logicalspike , [0; PS]', size(logicalspike,1), size(logicalspike,2) );
    lcif = lcif + lcif_ps;
end

%{
% NBCoupling 06-23-2014
if GLMType.CouplingFilters
    for pair=1:12
        lcif_cp = fastconvAH(pairspike{pair} , [0; CP{pair}]', size(pairspike{pair},1), size(pairspike{pair},2) );
        lcif = lcif + lcif_cp;
    end
end
% end NBCoupling
%}

%% Calculate some statistics about goodness of fit
if ~predict_only
    glm_ratepersec  = exp(lcif);
    glm_rateperbin  = params.bindur * glm_ratepersec;
    
    spikerate_bin    = size(find(logicalspike(:))) /  size(logicalspike(:));
    model_null0      = spikerate_bin * ones(1, params.bins);
    model_uop0       = (1/params.trials) * sum(logicalspike,1);
    model_null       = repmat(model_null0, params.trials, 1);
    model_uop        = repmat(model_uop0, params.trials, 1);
    null_logprob     = sum(eval_rasterlogprob(logicalspike, model_null, 'binary', 'conditioned'));
    uop_logprob      = sum(eval_rasterlogprob(logicalspike, model_uop, 'binary', 'conditioned'));
    % Check computations are correct %
    %null_logprob    = sum(eval_rasterlogprob(logicalspike, model_null0, 'notbinary', 'unconditioned'));
    %uop_logprob     = sum(eval_rasterlogprob(logicalspike, model_uop0, 'notbinary', 'unconditioned'));
    uop_bits          = uop_logprob - null_logprob;
    uop_bits_perspike = uop_bits / (sum(model_null0));
    uop_bits_persecond   = uop_bits / params.testdur_seconds;
    [raster_logprob_bin] = eval_rasterlogprob( logicalspike, glm_rateperbin,  'binary', 'conditioned') ;
    glm_logprob       = sum(raster_logprob_bin);
    glm_bits          = glm_logprob - null_logprob;
    glm_bits_perspike = glm_bits / (sum(model_null0));
    % glm_bits_perbin   = glm_bits / params.bins;
    glm_bits_persecond   = glm_bits / params.testdur_seconds;
    
    
    xvalperformance.logprob_null_raw     = null_logprob;
    xvalperformance.logprob_uop_raw      =  uop_logprob;
    xvalperformance.logprob_glm_raw      =  glm_logprob;
    xvalperformance.logprob_uop_bpspike  =  uop_bits_perspike;
    xvalperformance.logprob_glm_bpspike  =  glm_bits_perspike;
    xvalperformance.logprob_uop_bpsec    =  uop_bits_persecond;
    xvalperformance.logprob_glm_bpsec    =  glm_bits_persecond;
    xvalperformance.glm_normedbits       =  glm_bits_persecond / uop_bits_persecond;
    
end
%% Make the prediction rasters

lcif_const  = lcif_kx0 + lcif_mu0;
logical_sim = zeros(params.trials, params.bins);

if fittedGLM.GLMType.PostSpikeFilter && ~fittedGLM.GLMType.CouplingFilters
    cif_psgain = exp(PS);
    ps_bins     = length(cif_psgain);
    for i_trial = 1 : params.trials
        cif0         = exp(lcif_const);
        cif_ps       = cif0;
        binary_simulation = zeros(1,params.bins);
        for i = 1 : params.bins- ps_bins;
            roll = rand(1);
            if roll >  exp(-params.bindur*cif_ps(i));
                cif_ps(i+1: i + ps_bins) =  cif_ps(i+1: i + ps_bins) .* (cif_psgain');
                binary_simulation(i)= 1;
            end
        end
        logical_sim(i_trial,:) = binary_simulation ;
    end
    
    % NBCoupling
    %{
elseif GLMType.CouplingFilters && ~GLMType.PostSpikeFilter
    for i=1:12
        cif_cpgain{i} = exp(CP{i});
    end
    cp_bins     = length(cif_cpgain{1});
    for i_trial = 1 : size(pairspike{pair},1)
        cif0         = exp(lcif_const);
        cif_cp       = cif0;
        binary_simulation = zeros(1,params.bins);
        for i = 1 : params.bins- cp_bins;
            roll = rand(1);
            if roll >  exp(-params.bindur*cif_cp(i));
                binary_simulation(i)= 1;
            end
            for pair=1:12
                if pairspike{pair}(i_trial,i)
                    cif_cp(i+1: i + cp_bins) =  cif_cp(i+1: i + cp_bins) .* (cif_cpgain{pair}');
                end
            end
        end
        logical_sim(i_trial,:) = binary_simulation ;
    end
    
    
elseif GLMType.CouplingFilters && GLMType.PostSpikeFilter
    for i=1:12
        cif_cpgain{i} = exp(CP{i});
    end
    cif_psgain = exp(PS);
    ps_bins     = length(cif_psgain);
    cp_bins     = length(cif_cpgain{1});
    for i_trial = 1 : size(pairspike{pair},1)
        cif0         = exp(lcif_const);
        cif_ps_cp       = cif0;
        binary_simulation = zeros(1,params.bins);
        for i = 1 : params.bins- cp_bins;
            roll = rand(1);
            if roll >  exp(-params.bindur*cif_ps_cp(i));
                cif_ps_cp(i+1: i + ps_bins) =  cif_ps_cp(i+1: i + ps_bins) .* (cif_psgain');
                binary_simulation(i)= 1;
            end
            for pair=1:12
                if pairspike{pair}(i_trial,i)
                    cif_ps_cp(i+1: i + cp_bins) =  cif_ps_cp(i+1: i + cp_bins) .* (cif_cpgain{pair}');
                end
            end
        end
        logical_sim(i_trial,:) = binary_simulation ;
    end
    % end NBCoupling
    %}
else
    for i_trial = 1 : params.trials
        cif         = exp(lcif_const);
        binary_simulation = zeros(1,params.bins);
        for i = 1 : params.bins;
            roll = rand(1);
            if roll >  exp(-params.bindur*cif(i));
                binary_simulation(i)= 1;
            end
        end
        logical_sim(i_trial,:) = binary_simulation ;
    end
end

%% Save everything in the xval structure

if ~predict_only
    xvalperformance.rasters.recorded = logicalspike;
end
xvalperformance.rasters.glm_sim  = logical_sim;
xvalperformance.rasters.bintime  = params.bindur;
xvalperformance.rate  = cif_ps;

end





