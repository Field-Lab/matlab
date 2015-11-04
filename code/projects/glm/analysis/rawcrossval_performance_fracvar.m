function [crossval_score] = rawcrossval_performance_fracvar(raster_binned_sim,raster_binned_rec,smoothbins)
%%% PURPOSE %%%
% Compute cross-validated Fraction of Variance Explained
% Ability of Simuluated Raster to imitate (explain) recorded raster
% Only 


%%% INPUTS  %%%


%%% OUTPUTS %%%
% Crossval_Raw

%%% CALLS %%%
% none


% AKHeitman 2014-12-08 
%%

% BASIC QUANTITIES
reps        = size(raster_binned_rec,1); 
bins        = size(raster_binned_rec,2); 


% SETUP STRUCTURE
crossval_score.sigma_bin           = smoothbins;
crossval_score.metric_raw          = zeros(length(smoothbins),1);
crossval_score.metric_name         = 'Fraction of Variance Explained';
crossval_score.metric_note         = 'Uses just average signals, No Normalization';
crossval_score.notes.fracvar_eqn   = 'FracVar(model) := 1 - [L2norm(model - smoothed_spiketrain) / var(smoothed_spiketrain)]';
crossval_score.timestamp           = datestr(clock);
crossval_score.mfile_name          = mfilename('fullpath');


binarysignal_sim = sum(raster_binned_sim , 1) / reps;
binarysignal_rec = sum(raster_binned_rec , 1) / reps;


%%
% LOOP THROUGH THE DIFFERENT TIME SCALES
for i_smoothbin = 1:length(smoothbins)
    
    % DEFINE CONVOLUTION VECTOR
    sig_val        = smoothbins(i_smoothbin);
    convolve_index = [-4*sig_val:4*sig_val];
    convolve_vec   = normpdf(convolve_index,0,sig_val) / sum(normpdf(convolve_index,0,sig_val) );
    
    % SMOOTH EACH SIGNAL
    dummy_sim         = conv(binarysignal_sim, convolve_vec);
    dummy_rec         = conv(binarysignal_rec, convolve_vec);
    smoothsignal_sim  = dummy_sim((4*sig_val+1):(end-4*sig_val));
    smoothsignal_rec  = dummy_rec((4*sig_val+1):(end-4*sig_val));
    
    % FRAC VAR
    var_rec   = var(smoothsignal_rec);
    l2_error  = (1/bins) * sum( (smoothsignal_sim - smoothsignal_rec).^2);
    frac_var  = 1 - l2_error / var_rec;
    
    crossval_score.metric_raw(i_smoothbin) = frac_var;
end
    
end