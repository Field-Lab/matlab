function tc_fit_error = tc_fit_error(tc, fit_params, fixed_params, fit_indices, fixed_indices, verbose)
% tc_fit_error returns the error between a time course fit and a fit
%
% USAGE: tc_fit_error = tc_fit_error(tc, fit_params, fixed_params, fit_indices, fixed_indices, verbose)
%
% INPUTS:
%       sta             a spatial-timporal-chromatic STA
%       fit_params      parameter that will be fit
%       fixed_params    parameters that will be fixed
%       fit_indices     indices to the parameters that are fit
%       fixed_indices   indices to the parameters that are fixed
%       verbose         T or F that flags whether to plot spatial and temporal
%                       fits with the STA
%
% OPTIONAL INPUTS (none):
%       
%
% OUTPUTS:
%       fit_error       The RMSE between the fit and the STA
%
%
% NOTES:
%   This function is intended to be used as an anonymous function in
%   fminsearch or some other similar optimization funtion.  In fminsearch
%   the value of this function (the RMSE between time course and fit) will be
%   minimized.  ONLY BW TIME COURSE FUNCTIONS ARE SUPPORTED
%
% Author: GDF
% Date: 2014-05


% combine the fixed and free (fit) parameters
all_params(fit_indices) = fit_params;
all_params(fixed_indices) = fixed_params;

% get sta fit
tc_fit = time_course_fit_function(all_params);

% plot stuff if verbose is true
if verbose
    % spatial fit
    figure(1)
    
    % temporal fit
    norm_factor = max(abs(reshape(tc_fit, 1, [])));
    tc_fit = tc_fit ./ norm_factor;
    plot(tc_fit, '--k')
    hold on

    norm_factor = max(abs(reshape(tc, 1, [])));
    tc = tc ./ norm_factor;
    plot(tc, 'k')
    hold off
    drawnow
end

%tc_fit_error = sqrt(mean(tc - tc_fit').^2);
tc_fit_error = sqrt(mean((tc - tc_fit').^2));

end






