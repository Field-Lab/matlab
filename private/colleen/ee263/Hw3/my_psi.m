function y = my_psi(x,s)
    y = nan(length(x),1);
    for i = 1:length(x)
        if x(i)<-1
            y(i) = (x(i)+1)/s-1;
        elseif x(i)<=1
            y(i) = x(i);
        else
            y(i) = (x(i)-1)/s-1;
        end
    end
end