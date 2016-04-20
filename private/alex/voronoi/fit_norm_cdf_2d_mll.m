function [mllparams_x, fval_x, mllparams_y, fval_y] = fit_norm_cdf_2d_mll(x1, x2, y, start_points_x, start_points_y)

global xdata1 xdata2 ydata fact_precomp

xdata1 = x1';
xdata2 = x2';
ydata = y;
fact_precomp = factorial(ydata);

opts.MaxFunEvals=500;
opts.Display = 'off';

[mllparams_x, fval_x] = fminsearch(@tt_power, start_points_x, opts);
[mllparams_y, fval_y] = fminsearch(@tt1_power, start_points_y, opts);

% [mllparams_x, fval_x] = fminsearch(@tt1, start_points_x, opts);
% [mllparams_y, fval_y] = fminsearch(@tt2, start_points_y, opts);
% 
% [mllparams_x, fval_x] = fminsearch(@tt_log, start_points_x, opts);
% [mllparams_y, fval_y] = fminsearch(@tt1_log, start_points_y, opts);


% POWER
function y = tt_power(p) % x shift
global ydata xdata1 xdata2 fact_precomp

sat   = p(1);
b = p(2);
x0 = p(3);
yshift = p(4);
a = p(5);
y_interim = real((a*xdata1 + xdata2 - x0).^b);
y_interim = sat*(y_interim) + yshift;
y_interim(y_interim<=0) = 1e-6;
y_interim = y_interim .^ ydata .*exp(-y_interim) ./ fact_precomp;
y = -sum(log(y_interim));
if isinf(y);y = 10000000;end
if ~isreal(y);y = 10000000;end



function y = tt1_power(p) % y shift
global ydata xdata1 xdata2 fact_precomp

sat   = p(1);
b = p(2);
x0 = p(3);
yshift = p(4);
a = p(5);
y_interim = real((a*xdata1 - x0).^b) + real((xdata2 - x0).^b);
y_interim = sat*(y_interim) + yshift;
y_interim(y_interim<=0) = 1e-6;
y_interim = y_interim .^ ydata .*exp(-y_interim) ./ fact_precomp;
y = -sum(log(y_interim));
if isinf(y);y = 10000000;end
if ~isreal(y);y = 10000000;end


% LOGISTIC
function y = tt_log(p) % x shift
global ydata xdata1 xdata2 fact_precomp

sat   = p(1);
k = p(2);
mu = p(3);
c = p(4);
a = p(5);
y_interim = ones(size(xdata1,1),1).*sat./(1+ exp(-k*(a*xdata1 + xdata2 - mu))) + c;
y_interim(y_interim<=0) = 1e-6;
y_interim = y_interim .^ ydata .*exp(-y_interim) ./ fact_precomp;
y = -sum(log(y_interim));
if isinf(y);y = 1000000;end


function y = tt1_log(p) % y shift
global ydata xdata1 xdata2 fact_precomp

sat   = p(1);
k = p(2);
mu = p(3);
c = p(4);
a = p(5);
y_interim = ones(size(xdata1,1),1).*sat./(1+ exp(-k*(a*xdata1 - mu))) + ones(size(xdata2,1),1).*sat./(1+ exp(-k*(xdata2 - mu))) + c;
y_interim(y_interim<=0) = 1e-6;
y_interim = y_interim .^ ydata .*exp(-y_interim) ./ fact_precomp;
y = -sum(log(y_interim));
if isinf(y);y = 1000000;end


function y = tt1(p) % x shift
global ydata xdata1 xdata2 fact_precomp

sat   = p(1);
sigma = p(2);
mu = p(3);
sh = p(4);
a = p(5);
y_interim = sat .* normcdf(a*xdata1 + xdata2, mu, sigma) + sh;
y_interim(y_interim<=0) = 1e-6;
y_interim = y_interim .^ ydata .*exp(-y_interim) ./ fact_precomp;
y = -sum(log(y_interim));
if isinf(y);y = 10000000;end
if a>10;y = 10000000;end


function y = tt2(p) % y shift
global ydata xdata1 xdata2 fact_precomp

sat   = p(1);
sigma = p(2);
mu = p(3);
sh = p(4);
a = p(5);
y_interim = sat .* (normcdf(a*xdata1, mu, sigma) + normcdf(xdata2, mu, sigma)) + sh;
y_interim(y_interim<=0) = 1e-6;
y_interim = y_interim .^ydata .*exp(-y_interim) ./ fact_precomp;
y = -sum(log(y_interim));
if isinf(y);y = 10000000;end
if a>10;y = 10000000;end


