% AKHEITMAN 2014-11-12
% Actual Computation of rasters

% Expanded will include .. pairwise spike train differences ( embedded vs.
% viktor)
% variance in predictability of the spike train 

function raster_scores = raster_metrics_expanded(logicalspike,metricparams)
% Input metricparams  and logicalspike binary spike matrix dim (reps,bins)

% Unpack Parameters and Setup Raster Scores
pars.t_bin          = metricparams.bindur;
pars.sigma_bin(1)   = metricparams.smoothbins;
pars.sigma_bin(2)   = metricparams.slowsmoothbins;
pars.sigma_bin(3)   = metricparams.slowestsmoothbins;
pars.pair_numbers   = metricparams.pair_numbers;
pars.bins           = size(logicalspike,2);
pars.reps           = size(logicalspike,1);
pars.seconds        = pars.t_bin *size(logicalspike,2);
pars.spikes_pertrial= sum(logicalspike(:)) / pars.reps;
pars.spikes_persec  = pars.spikes_pertrial / ( pars.t_bin * size(logicalspike,2) );

raster_scores = cell(3,1); i_timescale = 1;
%%
for i_timescale = 1:3
   %% Set Time scale .. and parameters set
    sigma_bin             = pars.sigma_bin(i_timescale);
    convolve_index        = [-4*sigma_bin:4*sigma_bin];
    convolve_vec          = normpdf(convolve_index,0,sigma_bin) / sum(normpdf(convolve_index,0,sigma_bin) );
    
    raster_scores{i_timescale}.param_bintime   = metricparams.bindur;
    raster_scores{i_timescale}.param_sigma_bin = sigma_bin;
    raster_scores{i_timescale}.param_note = 'bintime: seconds, sigma_bin := standard deviation of smoothing gaussian';
    
    %% Log Prob and Bits Per Spike
    flatmodel.rate           = sum(logicalspike(:)) / (pars.reps*size(logicalspike,2)) ; 
    flatmodel.rate_rep       = flatmodel.rate*ones(size(logicalspike));
    flatmodel.logprob        = sum(eval_rasterlogprob(logicalspike, flatmodel.rate_rep, 'binary', 'conditioned'));
    %{ 
    %%% OLD COMPUTATION %%%
    prob_fulluop_bin      = sum(logicalspike,1) / (size(logicalspike,1));
    convolved_raw         = conv(prob_fulluop_bin, convolve_vec);
    model_uop_smooth00    = convolved_raw( (4*sigma_bin+1):(end-4*sigma_bin));
    model_uop_smooth      = repmat(model_uop_smooth00, reps, 1);
    logprob_smoothedrate  = sum(eval_rasterlogprob(logicalspike, model_uop_smooth, 'binary', 'conditioned'));
    %}
    convolved_raw_pertrial = zeros(size(logicalspike));
    for i_rep = 1:pars.reps
        convolved_raw                    = conv(logicalspike(i_rep,:), convolve_vec);
        convolved_raw_pertrial(i_rep,:)  = convolved_raw( (4*sigma_bin+1):(end-4*sigma_bin));
    end
    uop.model           = sum(convolved_raw_pertrial, 1) / pars.reps;
    [ ~, logprobmat]    = eval_rasterlogprob(logicalspike, repmat(uop.model, pars.reps, 1), 'binary', 'conditioned');
    uop.logprob_pertrial= sum(logprobmat,2);
    uop.logprob_pertrialpersec = uop.logprob_pertrial / pars.seconds;
    uop.bits_pertrial   = uop.logprob_pertrial - flatmodel.logprob; 
    uop.bps_pertrial    = uop.bits_pertrial / pars.spikes_pertrial;
    
    raster_scores{i_timescale}.bps_max  = max(uop.bps_pertrial);
    raster_scores{i_timescale}.bps_mean = mean(uop.bps_pertrial);
    raster_scores{i_timescale}.bps_std  = std(uop.bps_pertrial);
    raster_scores{i_timescale}.logprobpersec_max  = max(uop.logprob_pertrialpersec);
    raster_scores{i_timescale}.logprobpersec_mean = mean(uop.logprob_pertrialpersec);
    raster_scores{i_timescale}.logprobpersec_std  =  std(uop.logprob_pertrialpersec);
    
    
    %% Temporal Variation
    signal = uop.model;
    raster_scores{i_timescale}.temporal_variance_raw      = var(signal);
    raster_scores{i_timescale}.temporal_variance_normed   = var(signal) / (mean(signal)^2);
    raster_scores{i_timescale}.temporal_variance_normed_note = 'variance raw is variance, variance_normed is variance of smoothed firing rate divided by avg rate squared'; 
    
    %% Paired Measures, L2 norm and Viktor spike for trial variation
    raster_scores{i_timescale}.pairscores_note    = sprintf('based off of %d random combination of trials' , pars.pair_numbers); 
    raster_scores{i_timescale}.viktorspike_note   = sprintf('Cost computed with cost var ( 1/ (2*std)) where std is time of smoothing var' , pars.pair_numbers); 
    l2_diff = zeros(pars.pair_numbers,1);
    vksp    = zeros(pars.pair_numbers,1);
    viktor_cost = 1 / ( pars.t_bin * 2*sigma_bin );
    for i_pair = 1:pars.pair_numbers
        pair_vals = randsample(pars.reps,2);
        ind1 = pair_vals(1);
        ind2 = pair_vals(2);
        
        rate_1 = convolved_raw_pertrial(ind1,:);
        rate_2 = convolved_raw_pertrial(ind2,:);
        
        l2_diff(i_pair) = (1/pars.bins) * sum((rate_1 - rate_2).^2);
        
        spt_1  = pars.t_bin * find(logicalspike(ind1,:));
        spt_2  = pars.t_bin * find(logicalspike(ind2,:));
        
        vksp(i_pair) =spkd(spt_1, spt_2, viktor_cost);
    end
    raster_scores{i_timescale}.viktordist = mean(vksp);
    raster_scores{i_timescale}.viktordist_perspike = mean(vksp) / pars.spikes_pertrial;
    raster_scores{i_timescale}.square_error = mean(l2_diff);
    raster_scores{i_timescale}.fractional_error = mean(l2_diff) / var(signal) ; 
    
end



end