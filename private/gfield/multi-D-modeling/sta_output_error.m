function fit_error = sta_output_error(fit_params, linear_generator, spike_rate)

% report fit error for for a half squaring nonlinear output

nl_output = half_square_sta_nl(fit_params(1), fit_params(2), fit_params(3), linear_generator);

fit_error = sqrt(mean((spike_rate - nl_output).^2));




