function [params1, resnorm1, params2, resnorm2]=fit_norm_cdf_2d_lsq(x1,x2,y)

% bin y

x = [x1; x2]';


fr = zeros(1,14);
tmp_x=x1+x2;
tmp = linspace(min(tmp_x), max(tmp_x), 15);

for i = 2:length(tmp)
    fr(i-1) = mean(y(tmp_x>tmp(i-1) & tmp_x<=tmp(i)));
end
% 
% figure
% plot(tmp(1:end-1),fr)
% hold on
% 
% p=params1;
% xx = -0.3:0.1:0.3;
% % sat   = p(1);
% % b = p(2);
% % x0 = p(3);
% % yshift = p(4);
% sat   = 5.26;
% b = 3;
% x0 = -0.375;
% yshift = 0.2;
% yy = sat*((xx - x0).^3) + yshift;
% hold on
% plot(xx,yy)


% 
% xx = -0.5:0.1:0.5;
% sat   = 20.8;
% k = 4.3;
% mu = 0.94;
% a = 1;
% yy = ones(size(xx,1),1).*sat./(1+ exp(-k*(xx - mu)));
% hold on
% plot(xx,yy)
% tt = min(max(x1), max(x2));
% tt1 = max(min(x1), min(x2));
% tmp = linspace(tt1, tt, 7);
% sat   = params1(1);
% sigma = params1(2);
% mu = params1(3);
% sh = params1(4);
% a = params1(5);
% y = sat .* normcdf(a*tmp + tmp, mu, sigma) + sh;
% [xtt1, xtt2] = meshgrid(tmp,tmp);

   
% inits = [max(fr)-min(fr)    std(x1+x2)    mean(x1+x2)   min(fr) 1];

% inits = [1.8208    0.1292    0.2641    0.2666 1];


opts.MaxFunEvals=1500;
opts.Display = 'off';
% 
lb = [0 0 0 -Inf 0];
ub = [20 20 5 Inf Inf];
% inits = [(max(fr)-min(fr))*2    std(x1+x2)    0.5   min(fr) 1];
inits = [1.8208    0.1292    0.2641    0.2666 1];
fitfunc = @tt;
[params1, resnorm1] = lsqcurvefit(fitfunc, inits, x, y, lb,ub,opts);

lb = [0 0 0 -Inf 0];
ub = [50 50 10 Inf Inf];
% inits = [(max(fr)-min(fr))*5    std(x1+x2)*2    params1(3)*1.5   -params1(4) 1];
inits = [1.8208    0.1292    0.2641    0.2666 1];
fitfunc = @tt1;
[params2, resnorm2] = lsqcurvefit(fitfunc, inits, x, y, lb,ub,opts);
% 
% % 
% lb = [0 0 0 -1 0.1];
% ub = [50 20 5 1 10];
% % inits = [20    5    1    1 0];
% inits = [3  7  0.3  0 1];
% fitfunc = @tt_log;
% [params1, resnorm1] = lsqcurvefit(fitfunc, inits, x, y, lb,ub,opts);
% % [params1, resnorm1] = lsqcurvefit(fitfunc, inits, x, y, [],[],opts);
% 
% lb = [0 0 0 -1 0.1];
% ub = [50 20 5 1 10];
% inits = [10 8 0.5 0 1];
% fitfunc = @tt1_log;
% [params2, resnorm2] = lsqcurvefit(fitfunc, inits, x, y, lb,ub,opts);
% 
% lb = [0 0.1 -10 0 0];
% ub = [50 7 1 2 Inf];
% inits = [6 3 -0.36 min(fr) 1];
% fitfunc = @tt_power;
% [params1, resnorm1] = lsqcurvefit(fitfunc, inits, x, y, lb,ub,opts);
% % 
% lb = [0 0.1 -10 -2 0];
% ub = [10 10 1 2 Inf];
% inits = [1 5 -1 0 1];
% fitfunc = @tt1_power;
% [params2, resnorm2] = lsqcurvefit(fitfunc, inits, x, y, lb, ub,opts);


% POWER
function y = tt_power(p,x) % x shift
sat   = p(1);
b = p(2);
x0 = p(3);
yshift = p(4);
a = p(5);
tmp = real((a*x(:,1) + x(:,2) - x0).^b);
y = sat*(tmp) + yshift;
if ~isreal(sat)
    y = zeros(size(x(:,1)));
end


function y = tt1_power(p,x) % y shift
sat   = p(1);
b = p(2);
x0 = p(3);
yshift = p(4);
a = p(5);
tmp = real((a*x(:,1) - x0).^b) + real((x(:,2) - x0).^b);
y = sat*(tmp) + yshift;
if ~isreal(sat)
    y = zeros(size(x(:,1)));
end


% LOGISTIC
function y = tt_log(p,x) % x shift
sat   = p(1);
k = p(2);
mu = p(3);
c = p(4);
a = p(5);
y = ones(size(x,1),1).*sat./(1 + exp(-k*(a*x(:,1) + x(:,2) - mu))) + c;


function y = tt1_log(p,x) % y shift
sat   = p(1);
k = p(2);
mu = p(3);
c = p(4);
a = p(5);
y = ones(size(x,1),1).*sat./(1+ exp(-k*(a*x(:,1) - mu))) + ones(size(x,1),1).*sat./(1+ exp(-k*(x(:,2) - mu))) + c;


% CDF
function y = tt(p,x) % x shift
sat   = p(1);
sigma = p(2);
mu = p(3);
sh = p(4);
a = p(5);
y = sat .* normcdf(a*x(:,1) + x(:,2), mu, sigma) + sh;


function y = tt1(p,x) % y shift
sat   = p(1);
sigma = p(2);
mu = p(3);
sh = p(4);
a = p(5);
y = sat .* (normcdf(a*x(:,1), mu, sigma) + normcdf(x(:,2), mu, sigma)) + sh;





