function transforms2d(varargin)
%TRANSFORMS2D  description of functions operating on transforms
%
%   By 'transform' we mean an affine transform. A planar affine transform
%   can be represented by a 3x3 matrix.
%
%   Example
%   % create a translation by the vector [10 20]:
%   T = [1 0 10;0 1 20;0 0 1];
%
%   See also:
%   translation, rotation, scaling
%   homothecy, lineSymmetry
%   transformPoint, transformVector
%   transformLine, transformEdge
%
%
% ------
% Author: David Legland
% e-mail: david.legland@nantes.inra.fr
% Created: 2008-10-13,    using Matlab 7.4.0.287 (R2007a)
% Copyright 2008 INRA - BIA PV Nantes - MIAJ Jouy-en-Josas.
% Licensed under the terms of the LGPL, see the file "license.txt"

help('transforms2d');