function y = cumGauss(w, t, x)

% y = phi(z) = 0.5*(1+erf(z/sqrt(2)))
% z = w(x/t - 1)

z = w*((x/t) - 1);
y = 0.5*(1+erf(z/sqrt(2)));


