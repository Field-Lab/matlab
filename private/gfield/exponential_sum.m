function [fit] = exponential_sum(coef, x)

fit = coef(5) + coef(1) .* exp(-coef(2) .* x) + coef(3) .* exp(-coef(4) .* x);
