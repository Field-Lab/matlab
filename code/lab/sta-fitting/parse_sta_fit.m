function full_fit = parse_sta_fit(fit_struct)
%
% This function simply reorganizes the information in fit structure to be
% used in sta_fit_function.  Fit structure is a named set of parameters,
% while in sta_fit_function, those parameters are anonymous.
%

% 1 center_point_x
% 2 center_point_y
% 3 sd_x
% 4 sd_y
% 5 rotation_angle
% 6 color_weight_a
% 7 color_weight_b
% 8 color weight_c
% 9 x_dim
% 10 y_dim
% 11 num_colors
% 12 fit_surround
% 13 surround_sd_scale
% 14 surround_amp_scale (relative to center_amp)
% 15 scale_one
% 16 scale_two
% 17 tau_one
% 18 tau_two
% 19 n-filters
% 20 frame_number


params(1) = fit_struct.center_point_x;
params(2) = fit_struct.center_point_y;
params(3) = fit_struct.center_sd_x;
params(4) = fit_struct.center_sd_y;
params(5) = fit_struct.center_rotation_angle;
params(6) = fit_struct.color_weight_a;
params(7) = fit_struct.color_weight_b;
params(8) = fit_struct.color_weight_c;
params(9) = fit_struct.x_dim;
params(10) = fit_struct.y_dim;
params(11) = fit_struct.num_colors;
params(12) = fit_struct.initial_params(12);
params(13) = fit_struct.surround_sd_scale;
params(14) = fit_struct.surround_amp_scale;
params(15) = fit_struct.scale_one;
params(16) = fit_struct.scale_two;
params(17) = fit_struct.tau_one;
params(18) = fit_struct.tau_two;
params(19) = fit_struct.n_filters;
params(20) = fit_struct.frame_number;

full_fit = sta_fit_function(params);




