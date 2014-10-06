function polynomialCurves2d(varargin)
%POLYNOMIALCURVES2D  description of functions operating on polynomial curves
%
%   Polynomial curves are plane curves whose points are defined by a
%   polynomial for each x and y coordinate.
%   A polynomial curve is represented by 3 row vectors:
%   - the bounds of the parametrization
%   - the coefficients for the x coordinate (in increasing degree)
%   - the coefficients for the y coordinate (in increasing degree)
%
%   C = {[0 1], [5 5], [0 1 -1]};
%
%   As each coordinate are given by polynoms, its is possible to compute
%   various parameters like curvature, normal, of the exact length.
%
%   See also
%   polynomialCurvePoint, polynomialCurvePosition
%   polynomialCurveDerivative, polynomialCurveNormal,
%   polynomialCurveCurvature, polynomialCurveCurvatures
%   polynomialCurveProjection
%   polynomialCurveLength, polynomialCurveCentroid
%   polynomialCurveFit, polynomialCurveSetFit, polyfit2
%   polynomialDerivate
%
%
% ------
% Author: David Legland
% e-mail: david.legland@nantes.inra.fr
% Created: 2008-10-13,    using Matlab 7.4.0.287 (R2007a)
% Copyright 2008 INRA - BIA PV Nantes - MIAJ Jouy-en-Josas.
% Licensed under the terms of the LGPL, see the file "license.txt"

help('polynomialCurves2d');
