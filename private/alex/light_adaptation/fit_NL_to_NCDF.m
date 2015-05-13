function [p resnorm residual] = fit_NL_to_NCDF(xdata,ydata)

n = size(xdata,1);

% Set fitting bounds
lb = zeros(1,n+2);
ub = Inf(1,n+2);

% initials
psat0 = max(ydata(:));
sigma0 = 0.5;
p0 = [zeros(1,size(ydata,1)) psat0 sigma0];

[p resnorm residual] = lsqcurvefit(@normcdffitn, p0, xdata, ydata);%, lb, ub, 'MaxFunEvals', 5000, 'MaxIter', 5000);


function y = normcdffitn(p,x)
n = size(x,1);
sat   = p(n+1);
sigma = p(n+2);
for i = 1:n
    xscale = p(i);
    y(i,:) = sat .* normcdf(abs(x(i,:)) * xscale, 1, sigma);
end