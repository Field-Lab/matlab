function lines2d(varargin)
%LINES2D  description of functions operating on planar lines
%
%   A line is defined by a point (its origin), and a vector (its
%   direction). The different parameters are bundled into a row vector:
%   LINE = [x0 y0 dx dy];
%
%   A line contains all points (x,y) such that:
%   x = x0 + t*dx
%   y = y0 + t*dy;
%   for all t.
%
%   See also:
%   points2d, vectors2d, edges2d, rays2d
%   createLine, cartesianLine, medianLine
%   orthogonalLine, parallelLine, bisector
%   lineAngle, linePosition, intersectLines, projPointOnLine
%   onLine, distancePointLine, isLeftOriented
%   clipLine, invertLine
%   transformLine, drawLine
%   lineFit
%
%
% ------
% Author: David Legland
% e-mail: david.legland@nantes.inra.fr
% Created: 2008-10-13,    using Matlab 7.4.0.287 (R2007a)
% Copyright 2008 INRA - BIA PV Nantes - MIAJ Jouy-en-Josas.
% Licensed under the terms of the LGPL, see the file "license.txt"

help('lines2d');