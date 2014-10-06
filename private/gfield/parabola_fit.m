function y_vals = parabola_fit(coefs, x_vals)

y_vals = coefs(1).*x_vals.^2 + coefs(2);


