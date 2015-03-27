% Function that designs an interpolator for a given bandwidth and sampling rate

% Arguments

% tau - sampling rate (1 sample / tau bins) 
% T - maximum frequency (and effective length of signal)
% order - order of the filter (controls the steepness of the falloff btw pass and stop band in frequency domain )
function f = design_interpolator(tau,T,order)

if (nargin < 3)
    order = 3;
end

% Use Butterworth filter
[b a] = butter(order,tau/T);
in = [1; zeros(T-1,1)];
f = filter(b,a,in);

