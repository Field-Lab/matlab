function y = my_deriv(x, s, c)
    y = nan(1,3);
    y(1) = (1/c)*(x(1)-s(1))/norm(x(1:2)-s);
    y(2) = (1/c)*(x(2)-s(2))/norm(x(1:2)-s);
    y(3) = 1;
end
