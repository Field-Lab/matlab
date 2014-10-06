function sta_fit = sta_fit_from_obvius_fit(obvius_fit)
% sta_fit_from_obvius_fit     convert vision fit to sta fit in standard form
%
% usage:  sta_fit = sta_fit_from_obvius_fit(datarun,obvius_fit)
%
% arguments:    obvius_fit - struct with fields X
%
% outputs:       sta_fit - sta fit in standard form (mean, sd, angle, surround_sd_scale, surround_scale)
%
%
% 2009-04  gauthier
%
% NOTE: currently only spatial fit is supported, NOT time or color
%




% BODY OF THE FUNCTION

% initalize output
sta_fit = struct;

if isempty(obvius_fit) || isempty(fieldnames(obvius_fit))
    return
end

% convert center point
sta_fit.mean = obvius_fit.mean + [1 1];

% copy sd and angle
sta_fit.sd = obvius_fit.sd;
sta_fit.angle = obvius_fit.angle;
%fprintf('!!! OBVIUS FIT ANGLES MAY BE INCORRECT !!!\n')

% surround
sta_fit.surround_sd_scale = obvius_fit.surround_sd_scale;
sta_fit.surround_scale = obvius_fit.surround_scale;

