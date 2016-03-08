function [ptchs,ptchGrp] = createPatches(x,y,c, offset,FaceAlpha)
%createPatches.m
% This file will create a bar plot with the option for changing the
% FaceAlpha property. It is meant to be able to recreate the functionality
% of bar plots in versions prior to 2014b. It will create the rectangular
% patches with a base centered at the locations in x with a bar width of
% 2*offset and a height of y.

% Ensure x and y are numeric vectors
validateattributes(x,{'numeric'},{'vector'});
validateattributes(y,{'numeric'},{'vector'});
%#TODO Allow use of vector c

% Check size(x) is same as size(y)
assert(all(size(x) == size(y)),'x and y must be same size');

% Default FaceAlpha = 1
if nargin < 5
    FaceAlpha = 1;
end
if FaceAlpha > 1 || FaceAlpha <= 0
    warning('FaceAlpha has been set to 1, valid range is (0,1]');
    FaceAlpha = 1;
end

ptchGrp = hggroup; % For grouping the patches for legends

for k = 1:length(x)
    leftX = x(k) - offset; % Left Boundary of x
    rightX = x(k) + offset; % Right Boundary of x
    patch([leftX rightX rightX leftX],...
        [0 0 y(k) y(k)],c,'FaceAlpha',FaceAlpha,'Parent',ptchGrp);
end

ptchs = ptchGrp.Children;


end