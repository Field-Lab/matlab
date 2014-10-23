function fun = filtered_response(spk_times, tau)
%FILTERED_RESPONSE returns a function that represents a gaussian convolved with the spike train
%   fun = FILTERED_RESPONSE(spk_times, tau)
%   returns a function handle that accepts time as a [1 x n] vector and
%   returns the value of the spike train convolved with a gaussian at that
%   time where:
%       spk_times -> a [num_spikes x 1] vector of spike times
%       tau -> the tuning parameter to use.  The width of the gaussians.
% spk_times are a [nx1] array
%
%   For examples on usage see motion_script.m
%
% Marvin Thielk 2013
% mthielk@salk.edu

fun = @(t) sum(exp((-(repmat(spk_times, size(t)) - repmat(t, size(spk_times))) .^ 2) ./ (2 .* tau .^2)), 1);
end