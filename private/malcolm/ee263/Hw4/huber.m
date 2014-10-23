function y = huber(d, delta)
    y = nan(length(d),1);
    for i = 1:length(d)
        if abs(d(i)) <= delta
            y(i) = 0.5*d(i)^2;
        else
            y(i) = delta*(abs(d(i))-0.5*delta);
        end
    end
end