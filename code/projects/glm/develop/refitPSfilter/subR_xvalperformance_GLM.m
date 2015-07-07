function NL_xvalperformance     = subR_xvalperformance_GLM( stimdrivenrate, logicalspike, t_bin,PS)
% AKHEITMAN 2015-06-24  it works!
% New version includes PS filter

params.bindur     = t_bin;
params.bins       = length(stimdrivenrate);
params.trials     = size(logicalspike,1);
params.testdur_seconds = params.bindur * params.bins ;   

% Set log-conditional as stim driven only
lcif_teststim = log(stimdrivenrate);
lcif = repmat(lcif_teststim, params.trials,1);  

lcif_ps = fastconv(logicalspike , [0; PS]', size(logicalspike,1), size(logicalspike,2) );    
lcif = lcif + lcif_ps;

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
% PS Filter 
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

% IF NO PS FILTER
%{
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
%}
NL_xvalperformance.rasters.note           = 'glmsim includes altered non-linearity';
NL_xvalperformance.rasters.recorded       = logicalspike;
NL_xvalperformance.rasters.glm_sim        = logical_sim;
NL_xvalperformance.rasters.bintime        = params.bindur;
end              