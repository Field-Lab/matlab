function [L2_Error] = eval_xvalperformance_L2error(rastRec,rastSim, t_bin, smoothbins)

% AKHeitman 2014-11-24
% rastRec is a binary binned matrix of spikes recorded
% rastSim is the corresponding simulation designed to immitate
% t_bin is the size of each bin in the raster in seconds
% t_cost is the duration which will count for a unit of Viktor Spike

if (size(rastRec,1) ~= size(rastSim,1)) || (size(rastRec,2) ~= size(rastSim,2))
    error('rasters need to be of the same size!')
end
reps           = size(rastRec,1);
spikespertrial = length(find(rastRec(:)) == 1) / reps ;  
secspertrial   = t_bin * size(rastRec,2);
sigma_bin      = smoothbins;


convolve_index  = [-4*sigma_bin:4*sigma_bin];
convolve_vec    = normpdf(convolve_index,0,sigma_bin) / sum(normpdf(convolve_index,0,sigma_bin) );


ratebinRec      = sum(rastRec,1);
ratebinSim      = sum(rastSim,1);

ratesmoothRec   = conv(ratebinRec, convolve_vec) / reps;
ratesmoothSim   = conv(ratebinSim, convolve_vec) / reps;
ratesmoothRec   = ratesmoothRec( (4*sigma_bin+1):(end-4*sigma_bin));
ratesmoothSim   = ratesmoothSim( (4*sigma_bin+1):(end-4*sigma_bin));


signal_variance    = var(ratesmoothRec);
l2_error           = (1 / size(rastRec,2)) * sum ( (ratesmoothRec - ratesmoothSim).^2 );

L2_Error.raw = l2_error;
L2_Error.permillisec = 1000*l2_error / secspertrial;
L2_Error.perspikepermillisec    = 1000*l2_error / spikespertrial;
L2_Error.fractional  =  l2_error / signal_variance; 
L2_Error.signalvariance = signal_variance;
L2_Error.fractional_note = 'l2_error divided by variance of the signal';
end