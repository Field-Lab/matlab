delta = 1;
T = [t ones(N, 1)];
alpha_new = T\x; % initial guess is OLS solution (weights=1)
tol = 1e-3; % tolerance for convergence
dif = 1;
counter = 0;
while tol < dif
    alpha_current=alpha_new;
    w = huber(T*alpha_current-x,delta);
    w = w./w.^2;
    W = diag(sqrt(w));
    alpha_new = (W*T)\(W*x);
    dif = norm(alpha_new-alpha_current);
    counter = counter+1;
end

ols = T\x;
alpha = alpha_current;
plot(t,x,'x',t,ols(1)*t+ols(2),'-',t,alpha(1)*t+alpha(2),'--');