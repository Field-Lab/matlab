function fit_coef = fit_Gaussian(coef,x_points)

fit_coef = coef(1) * exp(-((x_points/(sqrt(2)*coef(2))).^2));