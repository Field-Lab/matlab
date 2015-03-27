function [paired_distance] = raster_paireddistances_vector(raster_binned,metric_params)
%%% PURPOSE %%%
% Compute inter-repition distance in raster trials
% Compute distance using a vector norm, smooth spike trains by time scales
% Average distance over a number of pairs

%%% INPUTS  %%%
% Raster_Binned
%  [trials,bins] matrix of zeros and ones
%  rows are different trials (blocks of spikes)
%  bins is just time bins
%  zeros = no spike , one = spike
% Metric_params
%  .smoothbins: smooth spike trains into "firing rates" with gaussians
%                of time scales of smooth bins
%  .pairnumbers: number of random pairs used to compute avg paired distance

%%% OUTPUTS %%%
% Paired_distance

%%% DETAILS %%%
% Using L2 Norm, with smoothed spike trains embedded into vector space
% Normalize by temporal variance of each signal
% L2 Norm ( rate1-rate2) / sqrt( var(rate1) * var(rate2)


%%% MAJOR CALLS %%%
% none

% AKHeitman 2014-12-05
% AKHeitman 2014-12-08 clean up / comment




%%
% UNPACK PARAMETERS
smoothbins  = metric_params.smoothbins;
pairnumbers = metric_params.pairnumbers;

% INTIAILIZE THE SOLUTION STRUCTURE
clear paired_distance
paired_distance.fracerr             = zeros(length(smoothbins),1);
paired_distance.smoothing_bins      = smoothbins;
paired_distance.fracerr_note        = sprintf('L2 Error between smoothed spiketrain, divided by variance of each smoothed spike train %d pairs' , pairnumbers);
paired_distance.timestamp           = datestr(clock);
paired_distance.mfile_name          = mfilename('fullpath');

% LOCAL PARAMS
pars.bins           = size(raster_binned,2);
pars.reps           = size(raster_binned,1);

% LOOP THROUGH TIMESCALES
for i_timescale = 1:length(smoothbins)
    % SMOOTH SPIKE TRAINS
    sigma_bin             = smoothbins(i_timescale);
    convolve_index        = [-4*sigma_bin:4*sigma_bin];
    convolve_vec          = normpdf(convolve_index,0,sigma_bin) / sum(normpdf(convolve_index,0,sigma_bin) ); 
    convolved_raw_pertrial = zeros(size(raster_binned));
    for i_rep = 1:pars.reps
        convolved_raw                    = conv(raster_binned(i_rep,:), convolve_vec);
        convolved_raw_pertrial(i_rep,:)  = convolved_raw( (4*sigma_bin+1):(end-4*sigma_bin));
    end

    % RUN THROUGH PAIRS, COMPUTE FRAC ERR 
    l2_diff_normvar = zeros(pairnumbers,1);
    for i_pair = 1:pairnumbers 
        pair_vals = randsample(pars.reps,2);
        ind1 = pair_vals(1);
        ind2 = pair_vals(2);
        rate_1 = convolved_raw_pertrial(ind1,:);
        rate_2 = convolved_raw_pertrial(ind2,:); 
        
        var1 = var(rate_1);
        var2 = var(rate_2);
        avg_var = .5*(var1 + var2);
        
        l2_diff_normvar(i_pair) = (1/pars.bins) * sum((rate_1 - rate_2).^2) / avg_var;
    end
    paired_distance.fracerr(i_timescale)  = mean(l2_diff_normvar) ; 
    
end

end