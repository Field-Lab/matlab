x_hat = zeros(2,T);
x_hat(:,1) = C\y(:,1); % estimate x_hat(1)
J = rho*norm(y(1)-C*x_hat(:,1))^2; % initialize J at t=1
F = [eye(n) ; sqrt(rho) * C];

% calculate each value of x_hat using previous one
for t = 2:T
    x_hat(:,t) = F\[A*x_hat(:,t-1)+B*u(:,t-1) ;
        sqrt(rho) * y(:,t)];
    J = J + norm(x_hat(:,t)-(A*x_hat(:,t-1)+B*u(:,t-1)))^2 ...
        + rho*norm(y(t)-C*x_hat(:,t))^2;
end

% J =
% 
%    7.1747e+05