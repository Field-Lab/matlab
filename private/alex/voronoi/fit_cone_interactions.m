function [params params1 resp resp1 x1, resp2, x2] = fit_cone_interactions(x, cone1, cone2)

x = x';

% 
% figure
% subplot(2,3,4)
% plot(cone1)
% hold on
% plot(cone2)

xt=[0.1 0.1 0.1];
[p] = lsqcurvefit(@(c,x) normcdf(x,c(1),c(2))+c(3),xt, x', cone1);
% subplot(2,3,5)
% r = -0.1:0.05:0.15;
% plot(r, normcdf(r,p(1),p(2))+p(3), '--k', 'linewidth', 2)
% hold on
% plot(x, cone1, 'linewidth', 2)

xt=[0];
[k, ress] = lsqcurvefit(@(t,x) normcdf(x,p(1),p(2))+p(3)+t(1),xt, x', cone2);
params = [p k];
resp = cone2-k(1);
% % figure
% plot(x,(cone2-k(1)), 'linewidth', 2)
% % legend('fit','no conditioning', 'conditioned on cone2')
% title(['y-shift fitting, r=', num2str(ress)])

%%%%%%%%%

xt=[0.1 0.1 0.1];
[p] = lsqcurvefit(@(c,x) normcdf(x,c(1),c(2))+c(3),xt, x', cone1);
% subplot(2,3,6)
% r = -0.4:0.05:0.11;
% plot(r, normcdf(r,p(1),p(2))+p(3), '--k', 'linewidth', 2)
% hold on
% plot(x, cone1, 'linewidth', 2)

% xt=[0,1];
% [k, ress] = lsqcurvefit(@(t,x) normcdf(x*t(2)+t(1),p(1),p(2))+p(3),xt, x', cone2);
% params1 = [p k];
% resp1 = cone2;
% x1 = x*k(2)+k(1);



xt=[0,1];
[k, ress] = lsqcurvefit(@(t,x) normcdf(x+t(1),p(1),p(2))+p(3),xt, x', cone2);
params1 = [p k];
resp1 = cone2;
x1 = x+k(1);


xt=[0,0];
[k, ress] = lsqcurvefit(@(t,x) normcdf(x+t(1),p(1),p(2))+p(3)+t(2),xt, x', cone2);
params1 = [p k];
resp2 = cone2-k(2);
x2 = x+k(1);
% figure
% plot(x*k(2)+k(1),cone2, 'linewidth', 2)
% % legend('fit','no conditioning', 'conditioned on cone2', 'location', 'best')
% title(['x-shift fitting, r=', num2str(ress)])