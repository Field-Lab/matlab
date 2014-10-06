function ell = inertiaEllipse(pts)
%INERTIAELLIPSE  inertia ellipse of a set of points
%   ELL = inertiaEllipse(PTS)
%
%   Example
%   inertiaEllipse
%
%   See also
%   circles2d
%
% ------
% Author: David Legland
% e-mail: david.legland@nantes.inra.fr
% Created: 2008-02-21,    using Matlab 7.4.0.287 (R2007a)
% Copyright 2008 INRA - BIA PV Nantes - MIAJ Jouy-en-Josas.
% Licensed under the terms of the LGPL, see the file "license.txt"
% affichage des labels individus

%TODO: take into account the orientation of the ellipse
center = mean(pts);
sigma  = var(pts);

ell = [center sigma*1.96 0];
