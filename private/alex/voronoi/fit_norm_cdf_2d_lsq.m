function [params1, resnorm1, params2, resnorm2]=fit_norm_cdf_2d_lsq(x1,x2,y)

% bin y

x = [x1; x2]';


fr = zeros(1,14);
tmp_x=x1+x2;
tmp = linspace(min(tmp_x), max(tmp_x), 15);

for i = 2:length(tmp)
    fr(i-1) = mean(y(tmp_x>tmp(i-1) & tmp_x<=tmp(i)));
end


   
% inits = [max(fr)-min(fr)    std(x1+x2)    mean(x1+x2)   min(fr) 1];

% inits = [1.8208    0.1292    0.2641    0.2666 1];


opts.MaxFunEvals=1500;
opts.Display = 'off';

lb = [0 0 0 -Inf 0];
ub = [20 20 5 Inf Inf];
inits = [(max(fr)-min(fr))*2    std(x1+x2)    0.5   min(fr) 1];
fitfunc = @tt;
[params1, resnorm1] = lsqcurvefit(fitfunc, inits, x, y, lb,ub,opts);

lb = [0 0 0 -Inf 0];
ub = [50 50 10 Inf Inf];
inits = [(max(fr)-min(fr))*3    std(x1+x2)    0.25   0 1];
fitfunc = @tt1;
[params2, resnorm2] = lsqcurvefit(fitfunc, inits, x, y, lb,ub,opts);


function y = tt(p,x) % x shift of first x data
sat   = p(1);
sigma = p(2);
mu = p(3);
sh = p(4);
a = p(5);
y = sat .* normcdf(a*x(:,1) + x(:,2), mu, sigma) + sh;


function y = tt1(p,x) % x shift of first x data
sat   = p(1);
sigma = p(2);
mu = p(3);
sh = p(4);
a = p(5);
y = sat .* normcdf(a*x(:,1), mu, sigma) + sat .* normcdf(x(:,2), mu, sigma) + sh;





