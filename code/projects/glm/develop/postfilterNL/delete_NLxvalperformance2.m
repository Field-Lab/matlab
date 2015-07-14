function NL_xvalperformance     = subR_xvalperformance_LNonly(cif_teststim, logicalspike, t_bin)
% AKHEITMAN 2015-06-24  it works!
params.bindur     = t_bin;
params.bins       = length(cif_teststim);
params.trials     = size(logicalspike,1);
params.testdur_seconds = params.bindur * params.bins ;   

% Set log-conditional as stim driven only
lcif_teststim = log(cif_teststim);
lcif = repmat(lcif_teststim, params.trials,1);  


% FOR NOW WE IGNORE PS FILTER  LN ONLY!!
%{
if PostSpikeFilter
    PS = otherfilters_NEW.PS; 
    lcif_ps = fastconv(logicalspike , [0; PS]', size(logicalspike,1), size(logicalspike,2) );    
    lcif = lcif + lcif_ps;
end
%}

glm_ratepersec  = exp(lcif);
glm_rateperbin  = params.bindur * glm_ratepersec;

spikerate_bin    = size(find(logicalspike(:))) /  size(logicalspike(:));      
model_null0      = spikerate_bin * ones(1, params.bins);
model_null       = repmat(model_null0, params.trials, 1);
null_logprob     = sum(eval_rasterlogprob(logicalspike, model_null, 'binary', 'conditioned'));
[raster_logprob_bin] = eval_rasterlogprob( logicalspike, glm_rateperbin,  'binary', 'conditioned') ;
glm_logprob       = sum(raster_logprob_bin);
glm_bits          = glm_logprob - null_logprob;
glm_bits_perspike = glm_bits / (sum(model_null0));
glm_bits_perbin   = glm_bits / params.bins;
glm_bits_persecond   = glm_bits / params.testdur_seconds;

NL_xvalperformance.note = 'Scores include optimized Non-Linearity';
NL_xvalperformance.logprob_null_raw            = null_logprob;
NL_xvalperformance.logprob_glm_raw      =  glm_logprob;
NL_xvalperformance.logprob_glm_bpspike  =  glm_bits_perspike;
NL_xvalperformance.logprob_glm_bpsec    =  glm_bits_persecond;



lcif_const  = lcif(1,:);
logical_sim = zeros(params.trials, params.bins);

% FOR NOW WE IGNORE PS FILTER  LN ONLY!!
%{
if PostSpikeFilter
    cif_psgain = exp(PS);
    ps_bins     = length(cif_psgain);
    for i_trial = 1 : size(logicalspike,1)
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
else
%}
for i_trial = 1 : size(logicalspike,1)
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

NL_xvalperformance.rasters.recorded       = logicalspike;
NL_xvalperformance.rasters.glm_withNL     = logical_sim;
NL_xvalperformance.rasters.bintime        = params.bindur;
end  




function NL_xvalperformance     = old_subR_NLxvalperformance(fittedGLM,lcif_teststim_NEW)
bpf               = fittedGLM.bins_per_frame;
params.bindur     = fittedGLM.t_bin;
params.bins       = length(lcif_teststim_NEW);
params.trials     = size(fittedGLM.xvalperformance.rasters.recorded,1);
params.testdur_seconds = params.bindur * params.bins ;   

logicalspike   = fittedGLM.xvalperformance.rasters.recorded;
raster_GLM_OLD = fittedGLM.xvalperformance.rasters.glm_sim;

% Set log-conditional as stim driven only
lcif_kx = repmat(lcif_teststim_NEW, params.trials,1);  
lcif    = lcif_kx;


% FOR NOW WE IGNORE PS FILTER  LN ONLY!!
%{
if fittedGLM.GLMType.PostSpikeFilter
    PS = otherfilters_NEW.PS; 
    lcif_ps = fastconv(logicalspike , [0; PS]', size(logicalspike,1), size(logicalspike,2) );    
    lcif = lcif + lcif_ps;
end
%}

glm_ratepersec  = exp(lcif);
glm_rateperbin  = params.bindur * glm_ratepersec;

spikerate_bin    = size(find(logicalspike(:))) /  size(logicalspike(:));      
model_null0      = spikerate_bin * ones(1, params.bins);
model_null       = repmat(model_null0, params.trials, 1);
null_logprob     = sum(eval_rasterlogprob(logicalspike, model_null, 'binary', 'conditioned'));
[raster_logprob_bin] = eval_rasterlogprob( logicalspike, glm_rateperbin,  'binary', 'conditioned') ;
glm_logprob       = sum(raster_logprob_bin);
glm_bits          = glm_logprob - null_logprob;
glm_bits_perspike = glm_bits / (sum(model_null0));
glm_bits_perbin   = glm_bits / params.bins;
glm_bits_persecond   = glm_bits / params.testdur_seconds;

NL_xvalperformance.logprob_null_raw            = null_logprob;
NL_xvalperformance.logprob_glm_withNL_raw      =  glm_logprob;
NL_xvalperformance.logprob_glm_withNL_bpspike  =  glm_bits_perspike;
NL_xvalperformance.logprob_glm_withNL_bpsec    =  glm_bits_persecond;

NL_xvalperformance.logprob_glm_raw     = fittedGLM.xvalperformance.logprob_glm_raw;
NL_xvalperformance.logprob_glm_bpspike = fittedGLM.xvalperformance.logprob_glm_bpspike;
NL_xvalperformance.logprob_glm_bpsec   = fittedGLM.xvalperformance.logprob_glm_bpsec;

lcif_const  = lcif_kx(1,:);
logical_sim = zeros(params.trials, params.bins);
if fittedGLM.GLMType.PostSpikeFilter
    cif_psgain = exp(PS);
    ps_bins     = length(cif_psgain);
    for i_trial = 1 : size(logicalspike,1)
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
else
    for i_trial = 1 : size(logicalspike,1)
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

NL_xvalperformance.rasters.recorded       = logicalspike;
NL_xvalperformance.rasters.glm_withNL     = logical_sim;
NL_xvalperformance.rasters.glm_original   = raster_GLM_OLD;
NL_xvalperformance.rasters.bintime        = params.bindur;
end 