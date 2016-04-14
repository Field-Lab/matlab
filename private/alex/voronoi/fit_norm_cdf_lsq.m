function [params_x, params_y]=fit_norm_cdf_lsq(x1,y1,x2,y2)

   
inits = [1.8208    0.1292    0.2641    0.2666    0];
opts.MaxFunEvals=1500;
opts.Display = 'off';

fitfunc = @tt;
params_x = lsqcurvefit(fitfunc, inits, [x1 x2], [y1 y2], [],[],opts);
fitfunc = @tt1;
params_y = lsqcurvefit(fitfunc, inits, [x1 x2], [y1 y2], [],[],opts);

function y = tt(p,x) % x shift of first x data
sat   = p(1);
sigma = p(2);
mu = p(3);
sh = p(4);
xshift = p(5);
x1 = x(:,1);
y(:,1) = (sat .* normcdf(x1, mu, sigma) + sh)';
x2 = x(:,2) + xshift;
y(:,2) = (sat .* normcdf(x2, mu, sigma) + sh)';


function y = tt1(p,x) % y shift of first y data
sat   = p(1);
sigma = p(2);
mu = p(3);
sh = p(4);
x1 = x(:,1);
y(:,1) = (sat .* normcdf(x1, mu, sigma) + sh)';
x2 = x(:,2);
yshift = p(5);
y(:,2) = (sat .* normcdf(x2, mu, sigma) + sh + yshift)';
