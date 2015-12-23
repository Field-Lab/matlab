function [probabilities, bins] = autocorrelation(spike_times, bin_size, duration)
% autocorrelation     calculates the autocorrelation in a spike train given
%                       a bin size and duration
%
% usage:  probabilities = autocorrelation(spike_times, bin_size, duration)
%
% arguments:     spike_times - spike times (in seconds)
%            bin_size - size of bins (in seconds)
%            duration - duration overwhich to calculate the autocorrelation
%                       in seconds
% outputs:     probabilities - normalized probability density across all
%                               spikes
%
%
% 2011-02 GDF
%


% BODY OF THE FUNCTION
 
num_spikes = length(spike_times);

% calculate interspike intervals;
spike_intervals = spike_times(2:num_spikes) - spike_times(1:num_spikes-1);

%set-up bins
bins = 0:bin_size:duration;
num_bins = length(bins);

% initialize probability vector
probabilities = zeros(1,num_bins); 

% calculate the mean spikes rate for "normalization"
%mean_rate = length(spike_times) ./ recording_length;

% calculate spike probablity for each bin
for bn = 1:(num_bins-1)
    temp_indices = find(spike_intervals > bins(bn) & spike_intervals < bins(bn+1));
    probabilities(bn) = length(temp_indices) ./ num_spikes;
end

