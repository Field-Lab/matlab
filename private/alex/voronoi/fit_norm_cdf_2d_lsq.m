function [params1, resnorm1, params2, resnorm2]=fit_norm_cdf_2d_lsq(x1,x2,y)

   
inits = [1.8208    0.1292    0.2641    0.2666 1 1];
opts.MaxFunEvals=1500;
opts.Display = 'off';

x = [x1; x2]';
fitfunc = @tt;
[params1, resnorm1] = lsqcurvefit(fitfunc, inits, x, y, [],[],opts);

x = [x1; x2]';
fitfunc = @tt1;
[params2, resnorm2] = lsqcurvefit(fitfunc, inits, x, y, [],[],opts);


function y = tt(p,x) % x shift of first x data
sat   = p(1);
sigma = p(2);
mu = p(3);
sh = p(4);
a = p(5);
b = p(6);
y = sat .* normcdf(a*x(:,1) + b*x(:,2), mu, sigma) + sh;


function y = tt1(p,x) % x shift of first x data
sat   = p(1);
sigma = p(2);
mu = p(3);
sh = p(4);
a = p(5);
b = p(6);
y = sat .* normcdf(a*x(:,1), mu, sigma) + sat .* normcdf(b*x(:,2), mu, sigma) + sh;




