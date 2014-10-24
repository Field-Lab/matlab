function [frame, tscale] = calc_reversing_grating_frame_intensities_jf(spec)
% WARNING: unconfirmed that this is actually what the stimulus computer
% puts out.  It's definitely close, but not guaranteed to be identical yet.

% This should now handle any orientation, but the above
% warning still applies! Orientaiton should be specified in degrees
% (0 is vertical, 90 is horizontal)

width  = spec.x_end - spec.x_start;
height = spec.y_end - spec.y_start;

orirad = (spec.orientation/360)*2*pi;
yinc = [1:height]*sin(orirad);
xinc = [1:width]*cos(orirad);
[yi xi] = ndgrid(yinc,xinc);
frame = sin((yi+xi+spec.spatial_phase)/spec.spatial_period*2*pi);

t = 0:(spec.temporal_period-1);
p = t / spec.temporal_period * 2 * pi;
tscale = sin(p);
