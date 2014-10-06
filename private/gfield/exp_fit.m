function y_vals = exp_fit(coef, x_vals)

y_vals = coef(1) .* exp(coef(2) .* x_vals) + coef(3);


