% XCORR_SPARSE   Compute cross-correlation between discrete time events
%
% usage:  T = xcorr_sparase(T1, T2, min, max, delta);
%
% arguments:      T1, T2 - column vector of events (must be int32)
%               min, max - range of correlation window
%                  delta - step size in output histogram
%
% outputs:  T - count histogram of correlation window
%
% Cross-correlation between discrete time events. This code is far faster
% than running XCORR on a "sparsified" version of the T1, T2. This is
% extremely useful for calculating correlation-based quantities in spike
% trains.
%
% All arguments must be integers. This makes the computation faster and
% not subject to round-off error.
%
% See also: COMPUTE_ACF, COMPUTE_CCF
%
% Mex file.
%
% shlens 2006-07-17
%
