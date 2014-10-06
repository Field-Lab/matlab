function boxes2d(varargin)
%BOXES2D  conventions for using bounding boxes
%
%   A box is represented as a set of limits in each direction:
%   BOX = [XMIN XMAX YMIN YMAX].
%
%   Boxes are used as result of computation for bounding boxes, and to clip
%   infinite shapes.
%
%   See also
%   clipLine
%   clipPolygon
%   polygonBounds
%
%
% ------
% Author: David Legland
% e-mail: david.legland@nantes.inra.fr
% Created: 2008-10-13,    using Matlab 7.4.0.287 (R2007a)
% Copyright 2008 INRA - BIA PV Nantes - MIAJ Jouy-en-Josas.
% Licensed under the terms of the LGPL, see the file "license.txt"

help('boxes2d');