function ray = createRay(varargin)
%CREATERAY create a ray (half-line)
%
%   RAY = createRay(POINT, ANGLE)
%   RAY = createRay(X0, Y0, ANGLE)
%   POINT is a N*2 array giving starting point of the ray, and ANGLE is the
%   orientation of the ray.
%
%   Ray is represented in a parametric form : [x0 y0 dx dy]
%   x = x0 + t*dx
%   y = y0 + t*dy;
%   for all t>0
%
%   See also:
%   rays2d, createLine, points2d
%
% ------
% Author: David Legland
% e-mail: david.legland@nantes.inra.fr
% Created: 2007-10-18
% Copyright 2007 INRA - BIA PV Nantes - MIAJ Jouy-en-Josas.
% Licensed under the terms of the LGPL, see the file "license.txt"

if length(varargin)==2
    p = varargin{1};
    theta = varargin{2};
    ray = [p cos(theta) sin(theta)];   
    
elseif length(varargin)==3   
    x = varargin{1};
    y = varargin{2};
    theta = varargin{3};
    ray = [x y cos(theta) sin(theta)];   

else
    error('Wrong number of arguments in ''createRay'' ');
end
