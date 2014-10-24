function [framesin, framecos, tdrift] = calc_drifting_grating_frame_intensities_jf(spec)
% WARNING: unconfirmed that this is actually what the stimulus computer
% puts out.  It's definitely close, but not guaranteed to be identical yet.

% This should now handle variable orientations, but the above
% warning still applies! Orientaiton should be specified in degrees,
% clockwise from the x-axis.

width  = spec.x_end - spec.x_start;
height = spec.y_end - spec.y_start;

orirad = (spec.orientation/360)*2*pi;
yinc = [1:height]*sin(orirad);
xinc = [1:width]*cos(orirad);
[yi xi] = ndgrid(yinc,xinc);

% for drifting gratings, we can generate a single frame each for
% sine and cosine, and then get any phase by taking a weighted sum
framesin = sin((yi+xi)/spec.spatial_period*2*pi);
framecos = cos((yi+xi)/spec.spatial_period*2*pi);

% range of phases between 0 and 2*pi
% for desired phase ph, compute framesin*cos(ph) + framecos*sin(ph)
tdrift = linspace(0,2*pi,spec.temporal_period);
tdrift = tdrift(1:end-1);
tdrift = tdrift(:);