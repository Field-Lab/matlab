function st = poisspikes(t, r)
% POISSPIKES
% usage: st = poisspikes(t, r)
%

% Expected # of spikes
expspikes = t*r;

% How many spikes we need to generate to be reasonably sure we'll fill the
% desired time window
sigma = 5; % Sigmas above t to aim for
nsafe = round(expspikes + sqrt(expspikes)*sigma);

% Generate exponentially distributed ISIs and then accumulate them for
% spiketimes
isis = -log(rand(nsafe,1))/r;
st = cumsum(isis);

if st(end) < t
    % Very improbably, we failed to generate enough spikes, so try again
    st = poisspikes(t,r);
else
    % Crop to time window
    st = st(st < t);
end