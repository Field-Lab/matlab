function rms = compute_rms(x,rms_mlength,rms_slength,padval)

if (~exist('padval','var'))
    padval = 0;
end

% Estimate the running mean
mean_est = running_mean(x,rms_mlength,padval);
% Estimate the running variance
var_est = (x - mean_est).^2;
% Estimate the root mean squared.
rms = sqrt(running_mean(var_est,rms_slength,mean(var_est(:))));

% Hacks 
%rms(abs(x) < eps) =1 ;
%rms (abs(rms) < eps) =1;


1;