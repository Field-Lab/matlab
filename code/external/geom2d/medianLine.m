function line = medianLine(varargin)
%MEDIANLINE create a median line between two points
%
%   L = medianLine(P1, P2);
%   Create the median line of points P1 and P2, that is the line containing
%   all points located at equal distance of P1 and P2.
%
%   L = medianLine(PTS);
%   Creates the median line of 2 points, given as a 2*2 array. Array has
%   the form:
%   [ [ x1 y1 ] ; [ x2 y2 ] ]
%
%   L = medianLine(EDGE);
%   Creates the median of the edge. Edge is a 1*4 array, containing [x0 y0]
%   coordinates of first point, and [dx dy], the vector from first point to
%   second point.
%  
%   See also:
%   lines2d, createLine, orthogonalLine
%
%   ---------
%
%   author : David Legland 
%   INRA - TPV URPOI - BIA IMASTE
%   created the 31/10/2003.
%

nargs = length(varargin);
x0=0;
y0=0;
dx=0;
dy=0;

if nargs == 1
    tab = varargin{1};
    if size(tab, 2)==2
        % array of two points
        x0 = tab(1,1); y0 = tab(1,2);
        dx = tab(2,1)-x0; dy = tab(2,2)-y0;
    else
        % edge equation
        x0 = tab(1); y0 = tab(2);
        dx = tab(3); dy = tab(4);
    end
elseif nargs==2
    % two points
    p1 = varargin{1};
    p2 = varargin{2};
    x0 = p1(1); y0 = p1(2);
    dx = p2(1)-x0; dy = p2(2)-y0;
end
    
        
line = [x0+dx/2, y0+dy/2, -dy, dx];