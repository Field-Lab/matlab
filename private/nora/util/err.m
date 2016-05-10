function error = err(a,b)
    %error = corr(a,b);
    error = sum(((a-b).^2))/sum(a);
end