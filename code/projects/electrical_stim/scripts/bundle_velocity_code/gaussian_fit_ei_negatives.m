function [fit, resnorm, residual, exitflag] = ...
    gaussian_fit_ei_negatives(ei, positions, A0, sigma0)

% Minimal electrodes are good initialization points.
[~, min_electrodes] = min(ei, [], 1);
    
% We only fit to negative signal.
ei(ei > 0) = 0;
    
% Fit params.
x0 = [A0, sigma0; positions(min_electrodes, :)];
[fit, resnorm, residual, exitflag] = ...
    lsqcurvefit(@symgaussgen, x0, positions, ei);