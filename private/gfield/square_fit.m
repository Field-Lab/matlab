function y_vals = square_fit(coef, x)

y_vals = (coef(1) .* x - coef(2)).^2;