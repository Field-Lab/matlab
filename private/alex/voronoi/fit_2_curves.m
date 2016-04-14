function [params_ref, params_x, params_y, params_xy, resnorm_x, resnorm_y, resnorm_xy] = fit_2_curves(datx,dat, sigm)

global tt

tt= sqrt(sigm);

dat = dat./tt;

opts = optimset('Display','off');

fitfunc = @normcdffitn_one;
inits = [0.15 0.1 0 0];
[params_ref, resn] = lsqcurvefit(fitfunc, inits,datx(:,1)', dat(:,1)',[],[],opts);
inits = [0.4 0.05 .05 min(dat(:,1))];
[params_ref2, resn2] = lsqcurvefit(fitfunc, inits,datx(:,1)', dat(:,1)',[],[],opts);

if resn2<resn
    params_ref = params_ref2;
end
% 
% x = datx(:,1);
% sat   = params_ref(1);
% sigma = params_ref(2);
% mu = params_ref(3);
% y = sat .* normcdf(x, mu, sigma);
% figure
% hold on
% plot(x,y, '-*')
% plot(datx(:,1),(dat(:,1)-0.39).*tt(:,1))
% plot(datx(:,1),(dat(:,2)-0.63).*tt(:,2))
% plot(datx(:,1),(dat(:,3)-0.2).*tt(:,3))
% 
% x = datx(:,1);
% sat   = params_y(4);
% sigma = params_y(5);
% mu = params_y(6);
% y = sat .* normcdf(x, mu, sigma);
% figure
% hold on
% plot(x,y, '-*')
% plot(datx(:,1),(dat(:,1)-params_y(1)).*tt(:,1))
% plot(datx(:,1),(dat(:,2)-params_y(2)).*tt(:,2))
% plot(datx(:,1),(dat(:,3)-params_y(3)).*tt(:,3))
% 
% 
% 
% fitfunc = @normcdffitn_y;
% inits = [-0.39 -0.63 -0.2 params_ref(1:3)];
% [params_y, resnorm_y] = lsqcurvefit(fitfunc, inits,datx', dat',[],[],opts);
% x = datx(:,1);
% sat   = params_y(4);
% sigma = params_y(5);
% mu = params_y(6);
% y = sat .* normcdf(x, mu, sigma);
% sum([((y+params_y(1))./tt(:,1)-dat(:,1)).^2; ...
%     ((y+params_y(2))./tt(:,2)-dat(:,2)).^2; ...
%     ((y+params_y(3))./tt(:,3)-dat(:,3)).^2])
% 
% 
% x = datx(:,1);
% sat   = params_ref(1);
% sigma = params_ref(2);
% mu = params_ref(3);
% y = sat .* normcdf(x, mu, sigma);
% sum([((y+0.39)./tt(:,1)-dat(:,1)).^2; ...
%     ((y+0.63)./tt(:,2)-dat(:,2)).^2; ...
%     ((y+0.2)./tt(:,3)-dat(:,3)).^2])
% 
% fitfunc = @normcdffitn_xy;
% inits = [zeros(1,2*size(datx,2)) params_ref(1:2)];
% [params_xy, resnorm_xy] = lsqcurvefit(fitfunc, inits,datx', dat',[],[],opts);
% 


% figure
% x = datx(:,1);
% sat   = 0.4;
% sigma = 0.05;
% mu = 0.05;
% sh = 0.4;
% y = (sat .* normcdf(x, mu, sigma)+sh)./tt(:,1);
% hold on
% plot(x,y, '-*')
% plot(datx(:,1),dat(:,1))
% 
% sum((y-dat(:,1)).^2)




fitfunc = @normcdffitn_x_mu;
inits = [zeros(1,size(datx,2)) params_ref];
[params_x, resnorm_x] = lsqcurvefit(fitfunc, inits,datx', dat',[],[],opts);
part=params_x;
part([1:3,end-1])=part([1:3,end-1])-params_x(end-1);
part(end-1) = [];
inits = part;
fitfunc = @normcdffitn_x;
[params_x, resnorm_x] = lsqcurvefit(fitfunc, inits,datx', dat',[],[],opts);

% fitfunc = @normcdffitn_x;
% inits = [zeros(1,size(datx,2)) params_ref([1 2 4])];
% [params_x, resnorm_x] = lsqcurvefit(fitfunc, inits,datx', dat',[],[],opts);


% fitfunc = @normcdffitn_y;
% inits = [zeros(1,size(datx,2)) params_ref(1:3)];
% [params_y, resnorm_y] = lsqcurvefit(fitfunc, inits,datx', dat',[],[],opts);


fitfunc = @normcdffitn_y_sh;
inits = [zeros(1,size(datx,2)) params_ref];
[params_y, resnorm_y] = lsqcurvefit(fitfunc, inits,datx', dat',[],[],opts);
part=params_y;
part([1:3,end])=part([1:3,end])+params_y(end);
part(end) = [];
inits = part;
fitfunc = @normcdffitn_y;
[params_y, resnorm_y] = lsqcurvefit(fitfunc, inits,datx', dat',[],[],opts);


fitfunc = @normcdffitn_xy_mu_sh;
inits = [0 params_x(1) 0 params_x(2) 0 params_x(3) params_x(4:6) 0];
[params_xy, resnorm_xy] = lsqcurvefit(fitfunc, inits,datx', dat',[],[],opts);
inits = [params_y(1) 0 params_y(2) 0 params_y(3) 0 params_y(4:5) 0 params_y(6)];
[params_xy2, resnorm_xy2] = lsqcurvefit(fitfunc, inits,datx', dat',[],[],opts);
if resnorm_xy>resnorm_xy2
    params_xy = params_xy2;
end

part=params_xy;
part([1:2:6,end])=part([1:2:6,end])+params_xy(end);
part([2:2:6,end-1])=part([2:2:6,end-1])-params_xy(end-1);
part(end-1:end) = [];

fitfunc = @normcdffitn_xy;
inits = part;
[params_xy, resnorm_xy] = lsqcurvefit(fitfunc, inits,datx', dat',[],[],opts);



function y = normcdffitn_one(params_ref,x)
global tt
sat   = params_ref(1);
sigma = params_ref(2);
mu = params_ref(3);
sh = params_ref(4);
y = (sat .* normcdf(x, mu, sigma)+sh)./tt(:,1)';


function y = normcdffitn_x(params_x,x)
global tt
n = size(x,1);
sat   = params_x(n+1);
sigma = params_x(n+2);
mu = 0;
sh = params_x(n+3);
for i = 1:n
    xscale = params_x(i);
    y(i,:) = (sat .* normcdf(x(i,:) + xscale, mu, sigma) + sh)./tt(:,i)';
end


function y = normcdffitn_x_mu(params_x,x)
global tt
n = size(x,1);
sat   = params_x(n+1);
sigma = params_x(n+2);
mu = params_x(n+3);
sh = params_x(n+4);
for i = 1:n
    xscale = params_x(i);
    y(i,:) = (sat .* normcdf(x(i,:) + xscale, mu, sigma) + sh)./tt(:,i)';
end

function y = normcdffitn_y(params_y,x)
global tt
n = size(x,1);
sat   = params_y(n+1);
sigma = params_y(n+2);
mu = params_y(n+3);
for i = 1:n
    yshift = params_y(i);
    y(i,:) = (sat .* normcdf(x(i,:), mu, sigma) + yshift)./tt(:,i)';
end


function y = normcdffitn_y_sh(params_y,x)
global tt
n = size(x,1);
sat   = params_y(n+1);
sigma = params_y(n+2);
mu = params_y(n+3);
sh = params_y(n+4);
for i = 1:n
    yshift = params_y(i);
    y(i,:) = (sat .* normcdf(x(i,:), mu, sigma) + yshift+sh)./tt(:,i)';
end

function y = normcdffitn_xy(params_xy,x)
global tt
n = size(x,1)*2;
sat   = params_xy(n+1);
sigma = params_xy(n+2);
mu = 0;
for i = 1:2:n
    yshift = params_xy(i);
    xscale = params_xy(i+1);
    y((i+1)/2,:) = (sat .* normcdf(x((i+1)/2,:) + xscale, mu, sigma) + yshift)./tt(:,(i+1)/2)';
end

function y = normcdffitn_xy_mu_sh(params_xy,x)
global tt
n = size(x,1)*2;
sat   = params_xy(n+1);
sigma = params_xy(n+2);
mu = params_xy(n+3);
sh = params_xy(n+4);
for i = 1:2:n
    yshift = params_xy(i);
    xscale = params_xy(i+1);
    y((i+1)/2,:) = (sat .* normcdf(x((i+1)/2,:) + xscale, mu, sigma) + yshift+sh)./tt(:,(i+1)/2)';
end
