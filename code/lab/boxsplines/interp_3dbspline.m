function ZI = interp_3dbspline(cx, cy, coeffs, XI, YI, order, nthreads)

if ~all(size(XI) == size(YI))
    error('XI and YI should be same size!');
end
ZI = zeros(size(XI));


for i = 1:numel(coeffs)
    coeff = coeffs(i);
    if coeff == 0
        continue;
    end
    
    x = (XI - cx(i)) ./ 30;
    y = (YI - cy(i)) ./ 30 * sqrt(3) / 2;
    
    if order == 2

        % Specialized code exists for this
        if nargin > 6
            % Multithreaded; may or may not be available
            BS = coeff*boxspline2_parallel(x, y, nthreads);
        else
            % Single threaded
            BS = coeff*boxspline2(x, y);
        end
        
    else
        
        % General code
        BS = coeff*boxsplinen(x, y, order);

    end
    ZI = ZI + BS;
end
