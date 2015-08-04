function t_frames = time_imrefresh_from_ttls(ttl_times, expected_ttl_interval)
%TIME_IMREFRESH_FROM_TTLS Timing of the image changes as recorded on the
%TTLS.
% 
%  T_FRAMES = TIME_IMREFRESH_FROM_TTLS(TTL_TIMES) calculates 
%  the time of the image refreshes in samples from a vector of ttl times 
%  also specified in samples. The k-th image appeared at sample
%  T_FRAMES(k).
%
%  TIME_IMREFRESH_FROM_TTLS(..., EXPECTED_TTL_INTERVAL) lets you specify
%  the expected number of samples between two triggers. It is used by the
%  function to estimate the image refresh rate in samples.
%
%  All indices follow matlab conventions and are 1-based.
%
%  Author: Georges Goetz - ggoetz@stanford.edu

% TODO: add in the possibility of specifying the following values as
% optional parameters (frames_per_ttl, expected_ttl_interval, 
% samples_per_frame)

frames_per_ttl = 100;

% When frames are continuously being updated, the duration between two TTLs
% should be at most the following value in samples.
if nargin == 1
    expected_ttl_interval = 0.95*20000; 
end

% Estimate the display refresh rate in samples
dttl = diff(ttl_times);
samples_per_frame = mean(dttl(dttl<expected_ttl_interval))/frames_per_ttl;

% Now, estimate the time in samples at which frames appeared
t_frames = [];
for kk = 1:length(ttl_times)
    t_frames = [t_frames ...
        int32(round((1:100)*samples_per_frame) + ttl_times(kk))]; %#ok<AGROW>
end

end % samples_to_images_from_ttls