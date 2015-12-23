function weights = spatial_diff_of_Gaussians(x, center_point, surround_sd_scale, surround_scale, locations)

%center_point(1) = x(1);
%center_point(2) = x(2);
center_scale = x(1);
center_sd = [abs(x(2)), 0; 0, abs(x(3))];
%surround_scale = abs(x(7));
angle = x(4);

% tranlate locations to polar coordinates
[thetas, rhos] = cart2pol(locations(1,:), locations(2,:));

% rotate coordinate frame by 'angle'
rotated_thetas = thetas - angle;

% translate back to cartesian coordinates
rotated_locations = locations;
[rotated_locations(1,:), rotated_locations(2,:)] = pol2cart(rotated_thetas, rhos);

%center_sd = [center_sd_x, center_sd_xy; center_sd_xy, center_sd_y];

center_weights = center_scale .* mvnpdf(rotated_locations, center_point, center_sd);

surround_weights = (center_scale * surround_scale) .* mvnpdf(rotated_locations, center_point, center_sd * surround_sd_scale);

weights= center_weights - surround_weights;








