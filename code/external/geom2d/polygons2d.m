function polygons2d(varargin)
%POLYGONS2D  description of functions operating on polygons
%
%   A simple polygon is represented by a N*2 array, each row of the array
%   representing the coordinates of a vertex. There is no need to repeat
%   the first point at the end.
%
%   For some functions, the orientation of the polygon can be relevant: CCW
%   stands for 'Conter-Clockwise' (positive orientation), CW stands for
%   'Clockwise'.
%
%   For multiple-connected polygons, the different connected boundaries are
%   separated by a row [NaN NaN].
%
%   Example:
%   % Simple polygon:
%   P1 = [1 1;2 1;2 2;1 2];
%   % Multiple polygon:
%   P2 = [10 10;40 10;40 40;10 40;NaN NaN;20 20;20 30;30 30;30 20];
%
%   See also
%   polygonCentroid, polygonArea, polygonLength, polygonBounds
%   polygonNormalAngle, polygonContains
%   clipPolygon, clipPolygonHP
%   splitPolygons, readPolygon, rectAsPolygon
%   steinerPoint, steinerPolygon, supportFunction, convexification
%   medialAxisConvex, polygonExpand
%   drawPolygon, fillPolygon
%
%
% ------
% Author: David Legland
% e-mail: david.legland@nantes.inra.fr
% Created: 2008-10-13,    using Matlab 7.4.0.287 (R2007a)
% Copyright 2008 INRA - BIA PV Nantes - MIAJ Jouy-en-Josas.
% Licensed under the terms of the LGPL, see the file "license.txt"

help('polygons2d');