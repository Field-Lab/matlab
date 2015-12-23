% problem 2b
% initialize alpha
alpha0 = [(max(x)-min(x))/2; 4; 0; mean(x)];
r0 = alpha0(1)*cos(2*pi*alpha0(2)*t+alpha0(3))+alpha0(4)-x;
% create iteration variables
Dr = nan(N, 4);
tol = 1e-3;
counter = 0;
alpha_new = alpha0;
alpha = zeros(4,1);
while norm(alpha-alpha_new) > tol
    alpha = alpha_new;
    Dr(:,1) = cos(2*pi*alpha(2)*t+alpha(3));
    Dr(:,2) = -2*pi*alpha(1)*t.*sin(2*pi*alpha(2)*t+alpha(3));
    Dr(:,3) = -alpha(1)*sin(2*pi*alpha(2)*t+alpha(3));
    Dr(:,4) = ones(N,1);
    r = alpha(1)*cos(2*pi*alpha(2)*t+alpha(3))+alpha(4)-x;
    alpha_new = alpha - Dr\r;
    counter = counter+1;
end

% alpha =
% 
%    22.4384
%     3.7400
%     0.9057
%     8.0294

% sum(r0.^2)
% 
% ans =
% 
%    4.2142e+04
% 
% sum(r.^2)
% 
% ans =
% 
%    6.7280e+03

xhat0 = x+r0;
xhat = x+r;
z = sortrows([t x xhat0 xhat],1);
plot(z(:,1),z(:,2),'.k','MarkerSize',14)
hold on
plot(z(:,1),z(:,3),':k',z(:,1),z(:,4),'-k')
xlabel('t')
ylabel('x')
legend({'x','xhat0','xhat'},'location','northeastoutside')
title('Problem 2b')