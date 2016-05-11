function F = symgaussgen(x, xdata)
% usage: F = symgaussgen(x, xdata)
%
% x holds the fit parameters.  The first row holds the sigma and scaling
% factor for the Gaussian.  After removing the first row, x is a t-by-2
% matrix giving x/y peak position for each time point.  xdata is a p-by-2
% matrix giving positions for each sampling point.  Output is a p-by-t
% matrix giving the symmetrical Gaussian signal at each sampling point
% for each time point.

A = x(1,1);
sigma = x(1,2);

xdist = bsxfun(@minus, x(2:end,1)', xdata(:,1));
ydist = bsxfun(@minus, x(2:end,2)', xdata(:,2));
dist = sqrt(xdist.^2 + ydist.^2);

F = A * normpdf(dist, 0, sigma);
