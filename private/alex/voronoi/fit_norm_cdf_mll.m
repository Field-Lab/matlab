function [mllparams_x, fval_x, mllparams_y, fval_y] = fit_norm_cdf_mll(x1, y1, x2, y2, start_points_x, start_points_y)

global xdata1 xdata2 ydata_com fact_precomp

xdata1 = x1;
xdata2 = x2;
ydata_com = [y1; y2]';
fact_precomp = factorial(ydata_com);

opts.MaxFunEvals=1500;
opts.Display = 'off';
[mllparams_x, fval_x] = fminsearch(@tt1, start_points_x, opts);
[mllparams_y, fval_y] = fminsearch(@tt2, start_points_y, opts);


function y = tt1(p) % x shift
global ydata_com xdata1 xdata2 fact_precomp

sat   = p(1);
sigma = p(2);
mu = p(3);
sh = p(4);
xshift = p(5);
x_com = [xdata1; xdata2 + xshift];
y_interim = (sat .* normcdf(x_com, mu, sigma) + sh)';
y_interim(y_interim<=0) = 1e-6;
y_interim = y_interim .^ydata_com .*exp(-y_interim) ./ fact_precomp;
y = -sum(log(y_interim));
if isinf(y);y = 1000000;end


function y = tt2(p) % y shift
global ydata_com xdata1 xdata2 fact_precomp

sat   = p(1);
sigma = p(2);
mu = p(3);
sh = p(4);
yshift = p(5);
y1 = (sat .* normcdf(xdata1, mu, sigma) + sh)';
y2 = (sat .* normcdf(xdata2, mu, sigma) + sh + yshift)';
y_interim = [y1 y2];
y_interim(y_interim<=0) = 1e-6;
y_interim = y_interim .^ydata_com .*exp(-y_interim) ./ fact_precomp;
y = -sum(log(y_interim));
if isinf(y);y = 1000000;end

