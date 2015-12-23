n = size(S,2);
x0 = zeros(3,1); % initial guess: zero
tol = 1e-3;
Dv = nan(n, 3);

x_current = x0;
x_new = x_current;
delta = 1;
counter = 0;
while delta>tol
    x_current = x_new;
    for i = 1:n
        Dv(i,:) = my_deriv(x_current, S(:,i), c);
    end
    x_new = x_current - pinv(Dv) * my_v(x_current,S,t,c);
    delta = norm(x_new - x_current);
    counter = counter+1;
end