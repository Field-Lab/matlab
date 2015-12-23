function xC = convx(xA, xB)
%CONVX
%   xC = conx(xA, xB);
%   Automatically generate the new x axis for a convolution.

error(nargchk(2, 2, nargin));

dxA = diff(xA);
dxB = diff(xB);

if isempty(dxA) 
    xC = xB;
    return
end

if any(((dxA - dxA(1)) / eps) > 5) || any(((dxB - dxB(1)) / eps) > 5)
    warning('CONVX:nonUniformXIncrements', 'Non-uniform x increments.  This may be due to round-off error, but be prepared for strange results.');
end
if ((dxA(1) - dxB(1)) / eps) > 5
    warning('CONVX:nonMatchingXIncrements', 'Non-matching x increments.  These functions should probably not be convolved.');
end

incr = dxA(1);
xCmin = min(xA) + min(xB);
xCmax = incr * (length(xA) + length(xB) - 2) + xCmin;

xC = xCmin:incr:xCmax;