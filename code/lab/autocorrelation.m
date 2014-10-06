function [probabilities, bins] = autocorrelation(spike_times, bin_size, duration, recording_length)
% autocorrelation     calculates the autocorrelation in a spike train given
%                       a bin size and duration
%
% usage:  [probabilities,bins] = autocorrelation(spike_times, bin_size, duration, recording_length)
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
% 2014-01 sravi modified ACF to calculate distribution of times between any 2 spikes
%               made code more efficient by vector indice increment,summing, looking at spike times in the duration range


% BODY OF THE FUNCTION

num_spikes = length(spike_times);

% calculate interspike intervals between all spikes heuristically;
k = 1;
 
for a = 2:num_spikes
    ST = spike_times(spike_times <spike_times(a-1,1)+duration+bin_size & spike_times > spike_times(a-1)) - spike_times(a-1);
    spike_intervals_all(1,k:k+(length(ST)-1)) = ST;
    k = k+length(ST);
end  


%set-up bins
bins = 0:bin_size:duration;
num_bins = length(bins);

% initialize probability vector
probabilities = zeros(1,num_bins); 

% calculate spike probablity for each bin
for bn = 1:(num_bins-1)
    probabilities(bn) = sum(spike_intervals_all > bins(bn) & spike_intervals_all < bins(bn+1));
end
probabilities(num_bins) = length(find(spike_intervals_all > bins(num_bins) & spike_intervals_all < (bins(num_bins)+bin_size)));

% calculate normalization parameters if required
% mean_rate = length(spike_times) ./ recording_length;
% mean_acf = mean(spike_intervals_all);
% var_acf = var(spike_intervals_all);
% rms_acf = rms(spike_intervals_all);
% norm = num_spikes;
% norm = recording_length;
% norm = mean_rate;
% norm = mean_acf;
% norm = var_acf;
% norm = sqrt(var_acf);
% norm = mean_acf*mean_rate;
% norm = rms_acf*mean_rate;
% norm = 1;
% norm = num_spikes/2
% probabilities = probabilities /2; %still do not know exact normalization vision does

end
