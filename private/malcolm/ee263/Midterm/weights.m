function w = weights(x,y,a)
    m = length(x);
    chi = [x ones(m,1)];
    delta = 1e-5;
    w = (abs(chi*a-y)./abs(y)) ./ max((chi*a-y).^2, delta);
end