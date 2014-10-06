function fit_error = sta_fit_error(sta, fit_params, fixed_params, fit_indices, fixed_indices, verbose)
% sta_fit_error returns the error between an STA fit and an STA
%
% USAGE: fit_error = sta_fit_error(sta, fit_params, fixed_params, fit_indices,fixed_indices, verbose)
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
%   the value of this function (the RMSE between STA and fit) will be
%   minimized.  
%
%   This function contains two subfunction. plot_spatial_sd is called if
%   verbose is true for plotting the Gaussian ellipse of the RF fit.
%   apply_all_constraints, ensures that various values of the fit remain
%   sensible.
%
% Author: GDF
% Date: 2011-06-11


% combine the fixed and free (fit) parameters
all_params(fit_indices) = fit_params;
all_params(fixed_indices) = fixed_params;

% get sta fit
sta_fit = sta_fit_function(all_params);

% plot stuff if verbose is true
if verbose
    % spatial fit
    figure(1)
    subplot(2,1,1)
    hold off
    temp_rf = rf_from_sta(sta);
    imagesc(norm_image(temp_rf))
    hold on
    plot_spatial_sd(all_params)
    axis image
    drawnow
    
    % temporal fit
   
    temp_stix = significant_stixels(sta, 'time', 'max', 'select', 'max', 'thresh', 3.5, 'robust_std_method', 1);
    fit_tc = time_course_from_sta(sta_fit, temp_stix);
    norm_factor = max(abs(reshape(fit_tc, 1, [])));
%     plot(fit_tc)
%     size(norm_factor)
%     size(fit_tc)
    fit_tc = fit_tc ./ norm_factor;
    subplot(2,1,2)
    if size(sta_fit, 3) == 3
        plot(fit_tc(:,1), '--r')
        hold on
        plot(fit_tc(:,2), '--g')
        plot(fit_tc(:,3), '--b')
    elseif size(sta_fit, 3) == 1
        plot(fit_tc, '--k')
        hold on
    else 
        error('dimensions of sta color is not recognized')
    end
    real_stix = significant_stixels(sta, 'time', 'max', 'select', 'max', 'thresh', 3.5, 'robust_std_method', 1);
    tc = time_course_from_sta(sta, real_stix);
    norm_factor = max(abs(reshape(tc, 1, [])));
    tc = tc ./ norm_factor;
    if size(sta_fit, 3) == 3
        plot(tc(:,1), 'r')
        plot(tc(:,2), 'g')
        plot(tc(:,3), 'b')
    elseif size(sta_fit,3) == 1
        plot(tc, 'k')
        hold off
    else
        error('dimensions of sta color is not recognized')
    end
    drawnow
end

% compute fit error
%fit_error = sqrt(mean((reshape(sta,1,[]) - reshape(sta_fit,1,[])).^2));

fit_error = sqrt(mean((reshape(sta,1,[]) - reshape(sta_fit,1,[])).^2)) +...
            (norm(all_params(6:8)) -1).^2; % this last bit constrains the color vector to be norm 1.

% % apply constraints
constraint_error = apply_constraints(all_params);
fit_error = fit_error + constraint_error;


end


% function to plot spatial fit ellispe
function plot_spatial_sd(params)

    % center
    ctr = [params(1) params(2)];
    rad_sd = [params(4) params(3)];
    rot_angle = params(5);
    
    % get points of an ellipse with these parameters
    [X, Y] = drawEllipse([ctr rad_sd (-1*rot_angle)]);
    plot(X, Y, 'k')
    
    % surround
    if params(13) ~= 0
        rad_sd_sur = rad_sd * params(13);
        [X, Y] = drawEllipse([ctr rad_sd_sur (-1*rot_angle)]);
        plot(X, Y, 'r')
    end
end



% function to constrain parameters
function constraint_error = apply_constraints(all_params)

    constraint_error = 0;

    % ensure surround_radius >= center_radius
    if all_params(13) < 1; constraint_error = constraint_error + 1000; end
    
    if all_params(13) > 5; constraint_error = constraint_error + 1000; end
    
    % ensure surround_scale >= 0
    if all_params(14) < 0;  constraint_error = constraint_error + 1000; end
    
    % ensure surround_scale < 1
    if all_params(14) >= 1;  constraint_error = constraint_error + 1000; end

    % ensure center_radius_x >= 0
    if all_params(3) < 0; 
        all_params(3) = realmin; 
        constraint_error = constraint_error + 1000;
    end

    % ensure center_radius_y >= 0
    if all_params(4) < 0
        all_params(4) = realmin;
        constraint_error = constraint_error + 1000;
    end
    
end

% % function to constrain parameters
% function all_params = apply_constraints(all_params)
% 
%     % ensure surround_radius >= center_radius
%     if all_params(13) < 1; all_params(13) = 1; end
%     
%     % ensure surround_scale >= 0
%     if all_params(14) < 0; all_params(14) = 0; end
%     
%     % ensure surround_scale < 1
%     if all_params(14) >= 1; all_params(14) = 1; end
% 
%     % ensure center_radius >= 0
%     if all_params(3) < 0; all_params(3) = realmin; end
% 
%     % ensure center_radius >= 0
%     if all_params(4) < 0; all_params(4) = realmin; end
%     
% end




