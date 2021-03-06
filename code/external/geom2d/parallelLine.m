function res = parallelLine(line, point)
%PARALLELLINE create a line parallel to another one.
%
%   RES = parallelLine(LINE, POINT);
%   Returns the line with same direction vector than LINE and going through
%   the point given by POINT. 
%   LINE is given as [x0 y0 dx dy] and POINT is [xp yp].
%
%
%   RES = parallelLine(LINE, DIST);
%   Uses relative distance to specify position. The new line will be
%   located at distance DIST, counted positive in the left side of LINE and
%   negative in the right side.
%
%   See also:
%   lines2d, orthogonalLine
%
%   ---------
%
%   author : David Legland 
%   INRA - TPV URPOI - BIA IMASTE
%   created the 31/10/2003.
%

%   HISTORY
%   31/07/2005 : add usage of distance

if size(point, 1)==1
    % use a distance. Compute position of point locatad at disatnce DIST on
    % the line orthogonal to the first one.
    point = pointOnLine([line(:,1) line(:,2) -line(:,4) line(:,3)], point);
end

% normal case : compute line through a point with given direction
res = zeros(1, 4);
res(1:2) = point;
res(3:4) = line(3:4);
