function [params_one, resnorm_one, params_one2, resnorm_one2, params, resnorm,paramsy, resnormy]=fit_one_poly(x1,y1,x2,y2)

if length(x1)>length(x2)    
    x = [x1(1:length(x2)) x2];
    y = [y1(1:length(y2)) y2];
else
    x = [x1 x2(1:length(x1))];
    y = [y1 y2(1:length(x1))];
end

fitfunc = @t;
inits = [108.4795   23.3808    3.8921    3.1341    0.7705];
opts.MaxFunEvals=500;
[params_one, resnorm_one] = lsqcurvefit(fitfunc, inits, x1', y1', [],[],opts);
[params_one2, resnorm_one2] = lsqcurvefit(fitfunc, inits, x2', y2', [],[],opts);

fitfunc = @tt;
inits = [119.0130   19.7748    2.6081    0.7839    0.3199  0];
opts.MaxFunEvals=1500;
[params, resnorm] = lsqcurvefit(fitfunc, inits, x', y', [],[],opts);
fitfunc = @tt1;
[paramsy, resnormy] = lsqcurvefit(fitfunc, inits, x', y', [],[],opts);

function y = t(p,x) % just one set of data

y = p(1)*x.^4 + p(2)*x.^3 + p(3)*x.^2 + p(4)*x + p(5);


function y = tt(p,x) % x shift of first x data

xscale = p(6);
x_t = x(1,:);
y(1,:) =  p(1)*x_t.^4 + p(2)*x_t.^3 + p(3)*x_t.^2 + p(4)*x_t + p(5);
x_t = xscale + x(2,:);
y(2,:) =  p(1)*x_t.^4 + p(2)*x_t.^3 + p(3)*x_t.^2 + p(4)*x_t + p(5);

function y = tt1(p,x) % y shift of first y data

yscale = p(6);
x_t = x(1,:);
y(1,:) =  p(1)*x_t.^4 + p(2)*x_t.^3 + p(3)*x_t.^2 + p(4)*x_t + p(5);
x_t = x(2,:);
y(2,:) =  p(1)*x_t.^4 + p(2)*x_t.^3 + p(3)*x_t.^2 + p(4)*x_t + p(5) + yscale;


