function err = spatial_diff_of_Gaussians_err(x, center_point, surround_sd_scale, surround_scale,locations, weights)

fit_weights = spatial_diff_of_Gaussians(x, center_point, surround_sd_scale, surround_scale, locations);

diff_weights = weights - fit_weights;

sqr_diff = diff_weights.^2;

summed_sqr_diff = sum(sqr_diff);

err = sqrt(summed_sqr_diff);
