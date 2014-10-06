function [probabilities, bins] = interspikeinterval(spike_times, bin_size, duration, recording_length)
% interspikeinterval     calculates the interspikeinterval in a spike train given
%                       a bin size and duration
%
% usage:  [probabilities,bins] = interspikeinterval(spike_times, bin_size, duration, recording_length)
%
% arguments:     spike_times - spike times (in seconds)
%            bin_size - size of bins (in seconds)
%            duration - duration overwhich to calculate the autocorrelation
%                       in seconds
%            recording_length = total length of stimulus
% outputs:     probabilities - normalized probability density across all
%                               spikes
%
%
% 2011-02 GDF
% 2014-04  sravi modified function name and function call from autocorrelation to interspikeintervals


% BODY OF THE FUNCTION
 
num_spikes = length(spike_times);

% calculate interspike intervals;
spike_intervals = spike_times(2:num_spikes) - spike_times(1:num_spikes-1);

%set-up bins
bins = 0:bin_size:duration;
num_bins = length(bins);

% initialize probability vector
probabilities = zeros(1,num_bins); 

% calculate the mean spike rate, mean and std dev of isi for "normalization" if required
% mean_rate = length(spike_times) ./ recording_length;
% mean_isi = mean(spike_intervals);
% var_isi = var(spike_intervals);

% calculate spike probablity for each bin
for bn = 1:(num_bins-1)
    temp_indices = find(spike_intervals > bins(bn) & spike_intervals < bins(bn+1));
    probabilities(bn) = length(temp_indices) ./ num_spikes;
end

end
