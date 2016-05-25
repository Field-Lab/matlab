function intersect = PD_GetIntersection( v0, vp, rmax, dr, f )
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This function evaluates points in range [v0-rmax;v0+rmax]
    % along vp direction, using 'f' function (step size: dr)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    r = -rmax:dr:rmax;
    %r = [r; r];
    wykres = zeros(length(r),4);
    %vr = ones(2,length(r)).*r;
    wykres(:,1) = r;
    wykres(:,2) = (r*vp(1) + v0(1))';
    wykres(:,3) = (r*vp(2) + v0(2))';
    wykres(:,4) = f(wykres(:,2),wykres(:,3));
    intersect = wykres;
end