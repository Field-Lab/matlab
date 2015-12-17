function [width] = inter_spike_interval(spikes)

times= nan(length(spikes)-1,1);
for i = 1:length(spikes)-1
    times(i) = spikes(i+1) - spikes(i);
end

bins = 0:0.001:max(times);

N = histc(times, bins);
N = N/sum(N);

width = fwhm(bins,N);

% initialize probability vector
