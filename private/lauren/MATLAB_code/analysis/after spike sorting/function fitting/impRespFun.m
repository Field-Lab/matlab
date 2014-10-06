function error = impRespFun(params, n, data)

% arguments:
%   params: parameters a and b to be fit
%
%   data: 2-d array, with data(1,:) vector of x-values and data(2,:) vector
%   of y-values
%
% returns: distance (error) between vector of data values and vector of calculated values based on current
% params


t0 = params(1);
tau = abs(params(2)); %to prevent tau from becoming negative
%tau = params;

t_sc = (data(1,:) - t0)/tau; %scaled/shifted time values
p1 = exp(-n*(t_sc-1));
p2 = t_sc.^n;
p = p1.*p2;

%zero out values corresponding to t_sc < 0
p(t_sc<0) = 0;

%t = data(1,:) - t0; %shifted t values
%p = (t.*exp(-t/tau)).^n;

%zero out any values for t<0



%find alpha that gives lowest error: Ax=b --> x = A\b
A = p';
b = data(2,:)';

alpha = A\b;


% figure
% hold on
% plot(data(2,:), 'k.')
% plot(alpha*p)
% hold off
% 
% 
% close(gcf)

error = norm(alpha*p - data(2,:));
