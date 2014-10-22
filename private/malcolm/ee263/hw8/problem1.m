% problem 1

% define problem data
A = [-0.1 1.0; -1.0 0.1];
c = [1 -1];
t_obs = [0.4 1.2 2.3 3.8];
y_obs = [1 -1 -1 1];

% at extrapolated time points t_e, draw points consistent with measurements
consistent = nan(4,1);
for t_e = [0.7 1.8 3.7]
    figure();
    plot([0 0] , [-1 +1] , 'k' , [-1 +1] , [0 0] , 'k');
    hold on;
    for x1 = -1:0.1:1
        for x2 = -1:0.1:1
            x = [x1 ; x2];
            for i = 1:4
                % compute extrapolated x and y
                x_e = expm((t_e-t_obs(i))*A)*x;
                y_e = sign(c*x_e);
                consistent(i) = (y_e == y_obs(i));
            end
            if all(consistent)
                plot(x1, x2, 'bo');
            end
            if sign(c*x) == 1
                plot(x1, x2, 'g^');
            else
                plot(x1, x2, 'rx');
            end
        end
    end
    hold off;
    title(sprintf('t_e = %d',t_e));
end