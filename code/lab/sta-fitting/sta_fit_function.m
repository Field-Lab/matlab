function full_fit = sta_fit_function(params)
%
% full_fit = sta_fit_function(params)
%
% This low level function receives a list of params (the order of which is
% important) and makes a spatial-temporal-chromatic STA based on those
% parameters.
%
%
% PARAM INFORMATION
% each number is an index into the vector 'params' and the name next to the
% number indicates how that value is used.
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
%
% Author: GDF 
% Data: 2011-06-10
%


% BODY OF FUNCTION

% RF center
spatial_filter = make_Gaussian_two_d('center_point_x', params(1),...
                        'center_point_y', params(2),...
                        'sd_x', params(3),...
                        'sd_y', params(4),...
                        'amp_scale', 1,...
                        'rotation_angle', params(5),...
                        'x_dim', params(9),...
                        'y_dim', params(10));


% RF surround
if params(12)
    surround_Gaussian = make_Gaussian_two_d('center_point_x', params(1),...
                            'center_point_y', params(2),...
                            'sd_x', params(3) * params(13),...
                            'sd_y', params(4) * params(13),...
                            'amp_scale', params(14),...
                            'rotation_angle', params(6),...
                            'x_dim', params(9),...
                            'y_dim', params(10));

    spatial_filter = spatial_filter - surround_Gaussian;

    % force peak of spatial filter to be 1 ( only need to do this if
    % surround is included in fit
    temp_sf_peak = max(reshape(spatial_filter, 1, []));
    spatial_filter = spatial_filter ./ temp_sf_peak;

end



% time course
%%%%
t_points = (1:1:params(20))-1;
t_filter_one = params(15) .* (t_points ./ params(17)).^params(19) .* exp(-params(19)*((t_points ./ params(17)) - 1));
t_filter_two = params(16) .* (t_points ./ params(18)).^params(19) .* exp(-params(19)*((t_points ./ params(18)) - 1));
tc = t_filter_one + t_filter_two;
tc = tc(params(20):-1:1);
%tc = abs(params(15)) .* (tc ./ abs(ext(tc)));

% figure(51)
% plot(tc)
% drawnow
%params(7)

% multiply the TC by the spatial-color RF and reshape
spatial_filter = reshape(spatial_filter, [], 1);
spatial_temporal_filter = spatial_filter * tc;
spatial_temporal_filter = reshape(spatial_temporal_filter, [], 1);

% add color information to fit
if params(11) == 1  % for BW
    full_fit = zeros(params(10), params(9), params(20), params(11));
    full_fit(:,:,:,1) = reshape(spatial_temporal_filter, [params(10), params(9), params(20), 1]);
    full_fit = permute(full_fit, [1,2,4,3]);
elseif params(11) == 3 % for RGB 
     color_vector = [params(6), params(7), params(8)];
     full_fit = spatial_temporal_filter * color_vector;
     full_fit = reshape(full_fit, [params(10), params(9), params(20), params(11)]);
     full_fit = permute(full_fit, [1,2,4,3]);
else
    error('num_color value is not supported')
end
   
% FIT INFORMATION 
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














