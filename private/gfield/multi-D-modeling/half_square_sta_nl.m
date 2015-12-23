function nl_output = half_square_sta_nl(fit_params, linear_generator)

nl_output = (fit_params(3) * (linear_generator - fit_params(1)).^2) - fit_params(2);
nl_output(nl_output < x_offset) = 0;