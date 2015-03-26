% AKHeitman 2014-11-26
% Re-evaluation of Raster Precision with goal of definitively identifying
% optimal percisino for each raster 

function [raster_precision] = rasterprecision_fullXVAL_old(logicalspike,smoothbins)
%%

% Grab mean and median values!! 
% BPS and UOP 
%%
%{
clear; clc
smoothbins = [.5 1 2 5 10 15 25 50 100];

load('dbug_rasterprecision.mat')

i_smoothbin = 1; i_rep = 1;
%}
%%
raster_precision.sigma_bin = smoothbins;
raster_precision.bps = zeros(length(smoothbins),1);
raster_precision.note_bps = 'Bits Per Spike : =  log prob of raster given time varying firing rate - log prob of raster avg firing rate ';
raster_precision.note_l2error = 'Error between avg firing rate and individual trial smoothed "firing rate", normed by variance of avg rate, or trial rate';
raster_precision.l2error_normmodel = zeros(length(smoothbins),1);
raster_precision.l2error_normtrial = zeros(length(smoothbins),1);

reps = size(logicalspike,1); bins = size(logicalspike,2); avgspikes = sum(logicalspike(:)) / reps;
for i_smoothbin = 1:length(smoothbins)
   %%
    sig_val        = smoothbins(i_smoothbin);
    convolve_index = [-4*sig_val:4*sig_val];
    convolve_vec   = normpdf(convolve_index,0,sig_val) / sum(normpdf(convolve_index,0,sig_val) );
    

    indtrials_smooth = zeros( reps, bins );
    for i_rep = 1:reps
        dummy = conv(logicalspike(i_rep,:), convolve_vec);
        indtrials_smooth(i_rep,:) = dummy((4*sig_val+1):(end-4*sig_val));
    end
    
    uop_bin        = sum(logicalspike,1) / (reps);
    uop_smooth     = sum(indtrials_smooth,1) / (reps);
    
    %%
    bps_vec = zeros(1,reps);
    error_normmodel_vec = zeros(1,reps);
    error_normtrial_vec = zeros(1,reps);
    for i_rep = 1:reps
        spiketrain = logicalspike(i_rep,:);
        
        model0 = indtrials_smooth;
        model0(i_rep,:) = [];
        model = sum(model0 ,1 ) / (reps-1);
        nullmodel = (1/bins) * sum(model) * ones(size(model));
        model = (1/reps) * ( (reps-1) * model + nullmodel); % hack to get no zero probabilities
        
        logprob_model     = sum(eval_rasterlogprob(spiketrain, model, 'binary', 'conditioned'));
        logprob_nullmodel = sum(eval_rasterlogprob(spiketrain, nullmodel, 'binary', 'conditioned'));
        bps_model = (logprob_model - logprob_nullmodel) / avgspikes;
        
        spiketrain_smooth = indtrials_smooth(i_rep,:);
        model_var = var(model);
        trial_var = var(spiketrain_smooth);
        l2_error  = (1/bins) * sum( (model-spiketrain_smooth).^2);
        
        bps_vec(i_rep) = bps_model;
        error_normmodel(i_rep) = l2_error / model_var;
        error_normtrial(i_rep) = l2_error / trial_var;
    end
    
    raster_precision.bps(i_smoothbin) = mean(bps_vec);
    raster_precision.l2error_normmodel(i_smoothbin) = mean(error_normmodel);
    raster_precision.l2error_normtrial(i_smoothbin) = mean(error_normtrial);

end
%%
    
end
    
    
 


