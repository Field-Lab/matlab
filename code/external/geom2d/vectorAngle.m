function alpha = vectorAngle(v, varargin)
%VECTORANGLE compute angle of a vector with horizontal axis
%
%   A = vectorAngle(V);
%   Returns angle between Ox axis and vector direction, in Counter
%   clockwise orientation.
%   The result is normalised between 0 and 2*PI.
%
%   See also:
%   vectors2d, vecnorm
%
% ------
% Author: David Legland
% e-mail: david.legland@nantes.inra.fr
% Created: 2007-10-18
% Copyright 2007 INRA - BIA PV Nantes - MIAJ Jouy-en-Josas.
% Licensed under the terms of the LGPL, see the file "license.txt"

alpha = mod(atan2(v(:,2), v(:,1))+2*pi, 2*pi);
