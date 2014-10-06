function weights = two_d_dog_fit(locations, fit_params)


center = fit_params.center;
center_scale = fit_params.center_scale;
center_sd = fit_params.center_sd;
angle = fit_params.angle;
surround_sd_scale = fit_params.surround_sd_scale;
surround_scale = fit_params.surround_scale;

    
num_points = size(locations,1);

covar_mat = diag(center_sd);
surround_covar_mat = covar_mat .* 2;

rotation_mat = [cos(angle), -sin(angle); sin(angle), cos(angle)];

centers = repmat(center, num_points, 1);
locations = (locations - centers) * rotation_mat;

weights = zeros(num_points,1);
for cn = 1:num_points

    center_gauss = exp(((locations(cn,:) / covar_mat) * locations(cn,:)') .* -0.5);

    surround_gauss = surround_scale .* exp(-((locations(cn,:) / surround_covar_mat) * locations(cn,:)') .* 0.5);
    
    weights(cn) = center_gauss + surround_gauss;
%    pause
end

% plot(weights)
% drawnow
% pause

