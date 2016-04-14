function [mleres_one1, mleres_one2, mleres1, mleres2] = get_mle(xdata,yy, xdata1, yy1, res_one1, res_one2, res1, res2)

global ydata x ydata1 x1

% Poisson
if length(xdata1)>length(xdata)
    ydata=yy;
    x = xdata;
    ydata1=yy1(1:length(xdata));
    x1 = xdata1(1:length(xdata));
else
    ydata=yy(1:length(xdata1));
    x = xdata(1:length(xdata1));
    ydata1=yy1;
    x1 = xdata1;
end


% opts.MaxFunEvals=1500;
% mleres = fminsearch(@tt, [res.p1, res.p2, res.p3, res.p4, res.p5], opts);

opts.MaxFunEvals=1500;
mleres_one1 = fminsearch(@t, res_one1, opts);
mleres_one2 = fminsearch(@t1, res_one2, opts);
mleres1 = fminsearch(@tt1, res1, opts);
mleres2 = fminsearch(@tt2, res2, opts);



function y = t(p)
global ydata x
k = p(1)*x.^4 + p(2)*x.^3 + p(3)*x.^2 + p(4)*x + p(5);
k(k<=0) = 1e-6;
y_prelog = k .^ydata .*exp(-k) ./ factorial(ydata);
y = -sum(log(y_prelog));
if isinf(y);y = 100000;end

function y = t1(p)
global ydata1 x1
k = p(1)*x1.^4 + p(2)*x1.^3 + p(3)*x1.^2 + p(4)*x1 + p(5);
k(k<=0) = 1e-6;
y_prelog = k .^ydata1 .*exp(-k) ./ factorial(ydata1);
y = -sum(log(y_prelog));
if isinf(y);y = 100000;end


function y = tt1(p)
global ydata x ydata1 x1

xscale = p(6);
x_scaled = x1 + xscale;
x_com = [x; x_scaled];
ydata_com = [ydata; ydata1];
k = p(1)*x_com.^4 + p(2)*x_com.^3 + p(3)*x_com.^2 + p(4)*x_com + p(5);
k(k<=0) = 1e-6;
y_prelog = k .^ydata_com .*exp(-k) ./ factorial(ydata_com);
y = -sum(log(y_prelog));
if isinf(y);y = 100000;end

function y = tt2(p)
global ydata x ydata1 x1

yshift = p(6);
ydata_com = [ydata; ydata1];
k1 = p(1)*x.^4 + p(2)*x.^3 + p(3)*x.^2 + p(4)*x + p(5);
k2 = p(1)*x1.^4 + p(2)*x1.^3 + p(3)*x1.^2 + p(4)*x1 + p(5) + yshift;
k = [k1; k2];
k(k<=0) = 1e-6;
y_prelog = k .^ydata_com .*exp(-k) ./ factorial(ydata_com);
y = -sum(log(y_prelog));
if isinf(y);y = 100000;end
