function R = rmatrix2(angle)
% RMATRIX   Generate 2-D rotation matrix
%
% Usage:  R = rmatrix2(angle)
%
% Generates a 2-D rotation matrix for rotating coordinates in a
% cartesian coordinate system (x,y) about the origin.
%
% Arguments:
%     angle - DEGREES rotation COUNTER-CLOCKWISE about origin
%

% assuming in degrees (ccw)
angle =- angle / 360 * 2*pi;

% make rotation matrix
R = [cos(angle) sin(angle); -sin(angle) cos(angle)];

