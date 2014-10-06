function circles2d(varargin)
%CIRCLES2D   description of functions operating on circles
%
%   Circles are represented by their center and their radius:
%   C = [xc yc r];
%   One sometimes considers orientation of circle, by adding an extra
%   boolean value in 4-th position, with value TRUE for direct (i.e.
%   turning Counter-clockwise) circles.
%
%   Circle arcs are represented by their center, their radius, the starting
%   angle and the angle extent.
%   CA = [xc yc r theta0 dtheta];
%   
%   Ellipses are represented by their center, their 2 semi-axis length, and
%   their angle with Ox direction.
%   E = [xc yc A B theta];
%
%
%   See also:
%   createCircle, createDirectedCircle, enclosingCircle
%   inCircle, onCircle
%   circleAsPolygon, circleArcAsCurve, ellipseAsPolygon
%   drawCircle, drawCircleArc, drawEllipse, drawEllipseArc
%
% ------
% Author: David Legland
% e-mail: david.legland@nantes.inra.fr
% Created: 2008-10-13,    using Matlab 7.4.0.287 (R2007a)
% Copyright 2008 INRA - BIA PV Nantes - MIAJ Jouy-en-Josas.
% Licensed under the terms of the LGPL, see the file "license.txt"

hrlp('circles2d');