function [str, flt_rsp1, flt_rsp2, flt_rsp1_shifted, flt_rsp2_shifted, spks_1_shifted, spks_2_shifted] = motion_signal(velocity, spks_1, spks_2, dx, trigger, trial_length, tau, tol)
%MOTION_SIGNAL calculates the motion signal for a given velocity between two spike trains
%   str = MOTION_SIGNAL(velocity, spks_1, spks_2, dx, trigger, trial_length, tau)
%   calculates the motion signal strength for the given velocity where:
%       spks_1 and spks_2 -> vertical vectors that contain all the spikes
%                            times for the two neurons. 
%       dx -> the distance between the two neurons (receptive fields)
%       trigger -> the start time for the trial you are using
%       trial_length -> length of the trial
%       tau -> the tuning parameter to use.  The width of the gaussians in
%              the filtered response algorithm
%
%   str = MOTION_SIGNAL(..., tol)
%   same as above with:
%       tol -> tolerance to use for the integration algorithm in matlab.
%              Defaults to 1e-3
%
%   [str, flt_rsp1, flt_rsp2, flt_rsp1_shifted, flt_rsp2_shifted, ...
%           spks_1_shifted, spks_2_shifted] = MOTION_SIGNAL(...)
%   returns some info useful for plotting and analyzing where:
%       flt_rsp1 and flt_rsp2 -> are function handles of the filtered
%                                responses of the spike trains. (functions
%                                of time best described over the interval
%                                of 0 <= t <= trial_length)
%       flt_rsp1_shifted and flt_rsp2_shifted -> are the funciton handles
%                                of the filtered responses to the shifted
%                                spike trains.
%       spks_1_shifted and spks_2_shifted -> are the spike times of the
%                                            shifted spikes (in vertical
%                                            vectors)
%
%   For examples on usage see motion_script.m
%
% Marvin Thielk 2013
% mthielk@salk.edu

if nargin < 8
    AbsTol = 1e-4;
    RelTol = 1e-3;
else
    RelTol = tol;
    AbsTol = tol * .1;
end

% trial start @ t=0
spks_1 = spks_1 - trigger;
spks_2 = spks_2 - trigger;
% only consider spikes that occured in the trial
spks_1 = spks_1(spks_1 >= 0 & spks_1 <= trial_length);
spks_2 = spks_2(spks_2 >= 0 & spks_2 <= trial_length);

if abs(dx / velocity) > trial_length / 2 ||isempty(spks_1) || isempty(spks_2)
    str = 0;
    flt_rsp1 = @(t) 0;
    flt_rsp2 = @(t) 0;
    flt_rsp1_shifted = @(t) 0;
    flt_rsp2_shifted = @(t) 0;
    spks_1_shifted = [];
    spks_2_shifted = [];
    return
end

% circularly shift spikes by dt
spks_1_shifted = spks_1 - dx / velocity;
spks_1_shifted = sort(mod(spks_1_shifted, trial_length));
spks_2_shifted = spks_2 - dx / velocity;
spks_2_shifted = sort(mod(spks_2_shifted, trial_length));
% replicate spikes before and after trial to minimize artifacts of spikes
% shifting circularly across the border
spks_1_shifted = [spks_1_shifted(ceil(end/2):end) - trial_length; spks_1_shifted; spks_1_shifted(1:floor(end/2)) + trial_length];
spks_2_shifted = [spks_2_shifted(ceil(end/2):end) - trial_length; spks_2_shifted; spks_2_shifted(1:floor(end/2)) + trial_length];
% filter responses
flt_rsp1 = filtered_response(spks_1, tau);
flt_rsp2 = filtered_response(spks_2, tau);
flt_rsp1_shifted = filtered_response(spks_1_shifted, tau);
flt_rsp2_shifted = filtered_response(spks_2_shifted, tau);
% and return the integral
str = integral(@(t) flt_rsp1(t) .* flt_rsp2_shifted(t) - flt_rsp2(t) .* flt_rsp1_shifted(t), 0, trial_length, ...
    'AbsTol', AbsTol, 'RelTol', RelTol);
%str = 0;
end