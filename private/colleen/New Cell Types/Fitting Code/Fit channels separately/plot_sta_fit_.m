function plot_sta_fit_(sta, fit, sig_stixels)
fit_params = fit{1}.fit_params;
fixed_params = fit{1}.fixed_params;
fit_indices = fit{1}.fit_indices;
fixed_indices = fit{1}.fixed_indices;

% combine the fixed and free (fit) parameters
all_params(fit_indices) = fit_params;
all_params(fixed_indices) = fixed_params;

% get sta fit
sta_fit = sta_fit_function(all_params);
    % spatial fit
    figure
    subplot(2,1,1)
    hold off
    temp_rf = rf_from_sta(sta);
    imagesc(norm_image(temp_rf))
    hold on
    plot_spatial_sd(all_params)
    axis image
    drawnow

 
    tc_full = time_course_from_sta(sta, sig_stixels);

%     norm_factor = max(abs(reshape(tc, 1, [])));
%     tc = tc ./ norm_factor;
  subplot(2,1,2)

        plot(linspace(1,all_params(20),all_params(20)),tc_full(:,1), 'r')
        hold on 
        plot(linspace(1,all_params(20),all_params(20)),tc_full(:,2), 'g')
        plot(linspace(1,all_params(20),all_params(20)),tc_full(:,3), 'b')
    drawnow

%% GREEN
    fit_params = fit{3}.fit_params;
fixed_params = fit{3}.fixed_params;
fit_indices = fit{3}.fit_indices;
fixed_indices = fit{3}.fixed_indices;

% combine the fixed and free (fit) parameters
all_params(fit_indices) = fit_params;
all_params(fixed_indices) = fixed_params;

% get sta fit
sta_fit = sta_fit_function(all_params);
   
    fit_tc = time_course_from_sta(sta_fit, sig_stixels);
    scale = max(tc_full(:,2));
%     norm_factor = max(abs(reshape(fit_tc, 1, [])));
%    fit_tc = fit_tc ./ norm_factor;
    subplot(2,1,2)
hold on

        plot(linspace(1,all_params(20),all_params(20)), fit_tc*scale/max(fit_tc), '--g')
      
 



%% RED  % Plot fit from whole curve fit because fit to just red is too noisy
    fit_params = fit{1}.fit_params;
fixed_params = fit{1}.fixed_params;
fit_indices = fit{1}.fit_indices;
fixed_indices = fit{1}.fixed_indices;

% combine the fixed and free (fit) parameters
all_params(fit_indices) = fit_params;
all_params(fixed_indices) = fixed_params;

% get sta fit
sta_fit = sta_fit_function(all_params);
    scale = max(tc_full(:,1));
    fit_tc = time_course_from_sta(sta_fit, sig_stixels);
%     norm_factor = max(abs(reshape(fit_tc, 1, [])));
%    fit_tc = fit_tc ./ norm_factor;
    subplot(2,1,2)
hold on
        plot(linspace(1,all_params(20),all_params(20)), fit_tc(:,1)*scale/max(fit_tc(:,1)), '--r')



%% BLUE
    fit_params = fit{4}.fit_params;
fixed_params = fit{4}.fixed_params;
fit_indices = fit{4}.fit_indices;
fixed_indices = fit{4}.fixed_indices;

% combine the fixed and free (fit) parameters
all_params(fit_indices) = fit_params;
all_params(fixed_indices) = fixed_params;

% get sta fit
sta_fit = sta_fit_function(all_params);
    scale = max(tc_full(:,3));
    fit_tc = time_course_from_sta(sta_fit, sig_stixels);
%     norm_factor = max(abs(reshape(fit_tc, 1, [])));
%    fit_tc = fit_tc ./ norm_factor;
    subplot(2,1,2)
hold on

        plot(linspace(1,all_params(20),all_params(20)), fit_tc*scale/max(fit_tc), '--b')
      
        
        plot(linspace(1,all_params(20),all_params(20)), zeros(1,all_params(20)), 'k')
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