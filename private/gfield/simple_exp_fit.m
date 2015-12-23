function y_vals = simple_exp_fit(coef, x)

y_vals = exp(coef(1) .* x - coef(2));