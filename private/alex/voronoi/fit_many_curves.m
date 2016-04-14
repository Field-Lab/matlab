function [p, g, qual, qual1] = fit_many_curves(datx,dat)


fitfunc = @normcdffitn;
inits = [zeros(1,size(datx,2)) 0.5    0.5    0.6];
inits = [zeros(1,size(datx,2)) 0.3686    0.2492   1];
[p, resnorm, residual] = lsqcurvefit(fitfunc, inits,datx', dat');
% 
% x = datx';
% n = size(x,1);
% sat   = p(n+1);
% sigma = p(n+2);
% mu = p(n+3);
% figure
% hold on
% for i = 1:n
%     xscale = p(i);
%     y(i,:) = sat .* normcdf(x(i,:) + xscale, mu, sigma);
%     plot(x(i,:),y(i,:))
% end
% figure
% plot(x,y)

qual = sum(residual.^2,2);

% fitfunc = @normcdffitn_one;
% inits = p(end-2:end);
% g_t = lsqcurvefit(fitfunc, inits,datx(:,1)', dat(:,1)');

% fitfunc = @normcdffitn_y;
% inits = [zeros(1,size(datx,2)) g_t];
% lb = [-inf(1,size(datx,2)) g_t*0.99];
% ub = [inf(1,size(datx,2)) g_t*1.01];
% [g, resnorm1, residual1] = lsqcurvefit(fitfunc, inits,datx', dat', lb, ub);

fitfunc = @normcdffitn_y;
inits = [zeros(1,size(datx,2)) p(end-2:end)];
[g, resnorm, residual] = lsqcurvefit(fitfunc, inits,datx', dat');

qual1 = sum(residual.^2,2);

function y = normcdffitn(p,x)
n = size(x,1);
sat   = p(n+1);
sigma = p(n+2);
mu = p(n+3);
for i = 1:n
    xscale = p(i);
    y(i,:) = sat .* normcdf(x(i,:) + xscale, mu, sigma);
end


% function y = normcdffitn_one(g,x)
% sat   = g(1);
% sigma = g(2);
% mu = g(3);
% y = sat .* normcdf(x, mu, sigma);


function y = normcdffitn_y(g,x)
n = size(x,1);
sat   = g(n+1);
sigma = g(n+2);
mu = g(n+3);
for i = 1:n
    yshift = g(i);
    y(i,:) = sat .* normcdf(x(i,:), mu, sigma) + yshift;
end