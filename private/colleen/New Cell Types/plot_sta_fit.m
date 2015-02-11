function plot_sta_fit(sta, fit_params, fixed_params, fit_indices, fixed_indices)
% combine the fixed and free (fit) parameters
all_params(fit_indices) = fit_params;
all_params(fixed_indices) = fixed_params;

% get sta fit
sta_fit = sta_fit_function(all_params);
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