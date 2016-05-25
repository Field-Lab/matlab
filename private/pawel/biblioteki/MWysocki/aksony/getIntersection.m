function intersect = getIntersection( v0, vp, rmax, dr, f )
    r = -rmax:dr:rmax;
    wykres = zeros(length(r),4);
    for i = 1:length(r)
       vr = v0 + r(i)*vp;
       wykres(i,1) = r(i);
       wykres(i,2) = vr(1); % x
       wykres(i,3) = vr(2); % y
       wykres(i,4) = f(vr(1),vr(2)); % value
    end
    intersect = wykres;
end

