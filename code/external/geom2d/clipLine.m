function edge = clipLine(line, box)
%CLIPLINE clip a line with a box
%
%   EDGE = clipLine(LINE, BOX);
%   LINE is a straight line given as a 4 element row vector: [x0 y0 dx dy],
%   with (x0 y0) being a point of the line and (dx dy) a direction vector,
%   BOX is the clipping box, given by its extreme coordinates: 
%   [xmin xmax ymin ymax].
%   The result is given as an edge, defined by the coordinates of its 2
%   extreme points: [x1 y1 x2 y2].
%   If line does not intersect the box, [NaN NaN NaN NaN] is returned.
%   
%   Function works also if LINE is a Nx4 array, if BOX is a Nx4 array, or
%   if both LINE and BOX are Nx4 arrays. In these cases, EDGE is a Nx4
%   array.
%      
%   See also:
%   lines2d
%   boxes2d
%   edges2d
%   clipEdge
%
% ------
% Author: David Legland
% e-mail: david.legland@nantes.inra.fr
% Created: 2007-08-27,    using Matlab 7.4.0.287 (R2007a)
% Copyright 2007 INRA - BIA PV Nantes - MIAJ Jouy-en-Josas.
% Licensed under the terms of the LGPL, see the file "license.txt"

%   HISTORY

% adjust size of two input arguments
if size(line, 1)==1
    line = repmat(line, size(box, 1), 1);
elseif size(box, 1)==1
    box = repmat(box, size(line, 1), 1);
elseif size(line, 1) ~= size(box, 1)
    error('bad sizes for input');
end

% allocate memory
edge = zeros(size(line, 1), 4);

% main loop on lines
for i=1:size(line, 1)
    % extract limits of the box
    xmin = box(i, 1);    xmax = box(i, 2);
    ymin = box(i, 3);    ymax = box(i, 4);
    
	% intersection with axis : x=xmin
	px1 = intersectLineEdge(line(i,:), [xmin ymin xmax ymin]);
	px2 = intersectLineEdge(line(i,:), [xmax ymin xmax ymax]);
	py1 = intersectLineEdge(line(i,:), [xmax ymax xmin ymax]);
	py2 = intersectLineEdge(line(i,:), [xmin ymax xmin ymin]);
	
	% sort points along the x coordinate, and  draw a line between
	% the two in the middle
	points = sortrows([px1 ; px2 ; py1 ; py2], 1);
	if points(2,1)>=xmin-1e-14 && points(2,1)<=xmax+1e-14
        if isfinite(points(3,1))
            edge(i, 1:4) = [points(2,:) points(3,:)];
        else
            edge(i, 1:4) = [points(1,:) points(2,:)];
        end 
    else
        % case of a line outside the box
        edge(i, 1:4) = [NaN NaN NaN NaN];
	end
end

