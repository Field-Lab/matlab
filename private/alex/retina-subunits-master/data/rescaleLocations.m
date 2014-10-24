function [out scaleVals] = rescaleLocations(locations)

% get argumnets
nCones = size(locations,1);
xx = locations(:,1);
yy = locations(:,2);
out = zeros(size(locations));

% get minimum values, ignoring 0s
mnx = min(xx(xx~=0));
mny = min(yy(yy~=0));

% subtract off minimums
xx = xx-mnx;
yy = yy-mny;

% get maximum value
mx = max(max(xx(:)),max(yy(:)));

% rescale to lie between 0 and 1
xx = xx/mx;
yy = yy/mx;

% rescale to have density of ~1 cones per square unit area
targetDensity = 1;
s = sqrt(nCones/targetDensity);

out(:,1) = xx*s;
out(:,2) = yy*s;

% operation was to subtract mnx (or mny), divide both by mx,
% and multiply by s
scaleVals.mnx = mnx;
scaleVals.mny = mny;
scaleVals.mx = mx;
scaleVals.s = s;
