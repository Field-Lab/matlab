function [raster_selfprediction] = raster_selfprediction_crossvalidated(raster_binned,smoothbins)

%%% PURPOSE %%%
% Compute cross-validated ability of raster to predict own spike trains
% Compute for different smoothing time scales
% Ultimately hope to find optimal time scale
% Ability of Avg Firing Rate to predict individual spike trains
% Both in Bits Per Spike and Fractional Variance


%%% INPUTS  %%%
% Raster_Binned
%  [trials,bins] matrix of zeros and ones
%  rows are different trials (blocks of spikes)
%  bins is just time bins
%  zeros = no spike , one = spike
% Smoothbins
%  smooth binary spike trains into "firing rates" with gaussians
%   of time scales of smooth bins

%%% OUTPUTS %%%
% Raster_Precision

%%% CALLS %%%
% eval_rasterlogprob

% AKHeitman 2014-11-26
% AKHeitman 2014-12-04 clean up / comment

%% 

% BASIC QUANTITIES
reps        = size(raster_binned,1); 
bins        = size(raster_binned,2); 
avgspikes   =  sum(raster_binned(:))/reps;

% SETUP STRUCTURE
raster_selfprediction.sigma_bin           = smoothbins;
raster_selfprediction.bps_mean            = zeros(length(smoothbins),1);
raster_selfprediction.fracvar_mean        = zeros(length(smoothbins),1);
raster_selfprediction.bps_std             = zeros(length(smoothbins),1);
raster_selfprediction.fracvar_std         = zeros(length(smoothbins),1);
raster_selfprediction.notes.crossval      = 'Each Computation is fully cross validated, If N reps,  than N different folds of raster';
raster_selfprediction.notes.scores        = '0 is null performance, 1 is perfect, negative means negative info';
raster_selfprediction.notes.model         = 'Time Varying firing rate of raster (without test trial), spikes smoothed by gaussian with standard deviations in bins';
raster_selfprediction.notes.rates         = 'Each spike convolved with Gaussian (up to 4std), where standard deviation is given by number of bins in sigma_bin';
raster_selfprediction.notes.bps_def       = 'BPS:= Bits Per Spike of the Model';
raster_selfprediction.notes.bps_eqn       = 'BPS(model):= logprob(spiketrain given model) - logprob(spiketrain given flat firing rate)';
raster_selfprediction.notes.fracvar_def   = 'FracVar:= Fraction of Variance of Smoothed Spike Train explained by Model';
raster_selfprediction.notes.fracvar_eqn   = 'FracVar(model) := 1 - [L2norm(model - smoothed_spiketrain) / var(smoothed_spiketrain)]';
raster_selfprediction.timestamp           = datestr(clock);
raster_selfprediction.mfile_name          = mfilename('fullpath');




% LOOP THROUGH THE DIFFERENT TIME SCALES
for i_smoothbin = 1:length(smoothbins)
    
    % DEFINE CONVOLUTION VECTOR
    sig_val        = smoothbins(i_smoothbin);
    convolve_index = [-4*sig_val:4*sig_val];
    convolve_vec   = normpdf(convolve_index,0,sig_val) / sum(normpdf(convolve_index,0,sig_val) );
    
    % SMOOTH EACH TRIAL
    indtrials_smooth = zeros( reps, bins );
    for i_rep = 1:reps
        dummy = conv(raster_binned(i_rep,:), convolve_vec);
        indtrials_smooth(i_rep,:) = dummy((4*sig_val+1):(end-4*sig_val));
    end
    
    % INITIALIZE SOLUTION 
    bps_vec     = zeros(1,reps);
    fracvar_vec = zeros(1,reps);
    
    % LOOP THROUGH EACH REP , FULLY CROSS VALIDATED
    for i_rep = 1:reps
        
        % TRIAL TO PREDICT
        spiketrain = raster_binned(i_rep,:);
        
        % CROSSVALIDATED MODEL
        model0 = indtrials_smooth;
        model0(i_rep,:) = [];
        model = sum(model0 ,1 ) / (reps-1);
        
        % HACK: GUARANTEE NO ZERO PROBABILITIES
        nullmodel = (1/bins) * sum(model) * ones(size(model));
        model     = (1/reps) * ( (reps-1) * model + nullmodel); % hack to get no zero probabilities
        
        % BITS PER SPIKE ON THE DISCRETE SPIKE TRAIN
        logprob_model     = sum(eval_rasterlogprob(spiketrain, model, 'binary', 'conditioned'));
        logprob_nullmodel = sum(eval_rasterlogprob(spiketrain, nullmodel, 'binary', 'conditioned'));
        bps_model = (logprob_model - logprob_nullmodel) / avgspikes;
        bps_vec(i_rep) = bps_model;
        
        % FRAC VAR ON THE SMOOTHED INDIVIDUAL FIRING RATE
        spiketrain_smooth = indtrials_smooth(i_rep,:);
        trial_var = var(spiketrain_smooth);
        l2_error  = (1/bins) * sum( (model-spiketrain_smooth).^2);
        fracvar_vec(i_rep) = 1 - l2_error / trial_var;
    end
    
    % FILL IN SOLUTION STRUCTURE
    raster_selfprediction.bps_mean(i_smoothbin)     = mean(bps_vec);
    raster_selfprediction.fracvar_mean(i_smoothbin) = mean(fracvar_vec);
    raster_selfprediction.bps_std(i_smoothbin)      = std(bps_vec);
    raster_selfprediction.fracvar_std(i_smoothbin)  = std(fracvar_vec);
end
    
end
    
    
 


