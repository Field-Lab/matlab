function v = my_v(x, s, t, c)
    v = nan(length(t),1);
    for i = 1:length(t)
        v(i) = (1/c)*norm(s(:,i)-x(1:2))+x(3)-t(i);
    end
end