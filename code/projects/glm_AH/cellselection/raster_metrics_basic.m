% AKHEITMAN 2014-11-12
% Actual Computation of rasters


function raster_scores = raster_metrics_basic(logicalspike,metricparams)

sigma_bin = metricparams.smoothbins;
t_bin = metricparams.bindur;
raster_scores.param_bintime  = t_bin;
raster_scores.param_sigma_bin      = sigma_bin;
raster_scores.param_note = 'bintime: seconds, sigma_bin := standard deviation of smoothing gaussian';


reps = size(logicalspike,1);

flatrate = sum(logicalspike(:)) / (reps*size(logicalspike,2)) ; 

convolve_index     = [-4*sigma_bin:4*sigma_bin];
convolve_vec       = normpdf(convolve_index,0,sigma_bin) / sum(normpdf(convolve_index,0,sigma_bin) );
prob_fulluop_bin   = sum(logicalspike,1) / (size(logicalspike,1));
convolved_raw      = conv(prob_fulluop_bin, convolve_vec);
model_uop_smooth0  = convolved_raw( (4*sigma_bin+1):(end-4*sigma_bin));
model_uop_smooth   = repmat(model_uop_smooth0, reps, 1);

flatmodel          = flatrate*ones(size(logicalspike));

logprob_smoothedrate  = sum(eval_rasterlogprob(logicalspike, model_uop_smooth, 'binary', 'conditioned'));
logprob_flatmodel     = sum(eval_rasterlogprob(logicalspike, flatmodel, 'binary', 'conditioned'));



spikes.pertrial = sum(logicalspike(:)) / reps;
spikes.persec   = spikes.pertrial / ( t_bin * size(logicalspike,2) );
raster_scores.spikespersec = spikes.persec;
% trial variability
spikecount_bytrial = sum(logicalspike,2);
raster_scores.trial_variance      = std(spikecount_bytrial) / mean(spikecount_bytrial)  ;
raster_scores.trial_variance_note = 'standard dev of spikes per trial devided by avergae spikes per trial'; 
% smoothed temporal variation
smoothed_rate = model_uop_smooth0;


raster_scores.temporal_variance      = std(smoothed_rate) / mean(smoothed_rate);
raster_scores.temporal_variance_note = 'standard dev of smoothed firing rate divided avg value of bin'; 
raster_scores.logprob.smoothedrate  = logprob_smoothedrate;
raster_scores.logprob.flatmodel     = logprob_flatmodel;
raster_scores.bits           = (logprob_smoothedrate - logprob_flatmodel);
raster_scores.bits_per_spike = (logprob_smoothedrate - logprob_flatmodel) / spikes.pertrial;
raster_scores.bits_per_sec   = (logprob_smoothedrate - logprob_flatmodel) / spikes.persec;

raster_scores.avgrate = smoothed_rate;




end