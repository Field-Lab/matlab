function points2d
%POINTS2D  description of functions operating on points
%
%   A point is defined by its two cartesian coordinate, put into a row
%   vector of 2 elements:
%   P = [x y];
%
%   Several points are stores in a matrix with two columns, one for the
%   x-coordinate, one for the y-coordinate.
%   PTS = [x1 y1 ; x2 y2 ; x3 y3];
%
%   Example
%   P = [5 6];
%
%   See also:
%   centroid, polarPoint
%   angle2Points, angle3Points, angleSort
%   distancePoints, minDistancePoints
%   pointOnLine
%   transformPoint, clipPoints, drawPoint
%
% ------
% Author: David Legland
% e-mail: david.legland@nantes.inra.fr
% Created: 2008-10-13,    using Matlab 7.4.0.287 (R2007a)
% Copyright 2008 INRA - BIA PV Nantes - MIAJ Jouy-en-Josas.
% Licensed under the terms of the LGPL, see the file "license.txt"

help('points2d');