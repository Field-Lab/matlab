function y = my_mu(x)
    y = nan(length(x),1);
    for i = 1:length(x)
        if x(i)<-1
            y(i) = -1;
        elseif x(i)<=1
            y(i) = x(i);
        else
            y(i) = 1;
        end
    end
end