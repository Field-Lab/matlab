function [mu, sigma, offset] = fitGaussianToRF(neuronRFdata)
% Fits a 2D gaussian to the eRF for a single neuron.

xPos = neuronRFdata.pixelPos(:,1);
yPos = neuronRFdata.pixelPos(:,2);
zPos = neuronRFdata.allResponses/max(neuronRFdata.allResponses);

% Starting point
x0 = sum(xPos.*zPos)/sum(zPos);
y0 = sum(yPos.*zPos)/sum(zPos);
a = 1e-10;
b = 0;
c = 1e-10;
d = 0;

% Setting options for the termination condition of the fit
fitoptions = optimset('TolFun',1e-6,'TolX',1e-5);

% Solving the non-linear least-squares problem
initPos = [a b c d x0 y0];
[xres, resnorm, residual, exitflag] = lsqcurvefit(@fitGauss2d, initPos, [xPos yPos], zPos, ...
    [0 0 0 0 min(xPos) min(yPos)], [100 100 100 100 max(xPos) max(yPos)], fitoptions);

mu = [xres(5) xres(6)];
sigma = inv([xres(1) xres(2); xres(2) xres(3)]);
offset = xres(4);

end % fitGaussianToERF

function val = fitGauss2d(x,XY)
% Function we're trying to fit to the data
% x: vector containing the different coefficients of the equation.
% XY = [X, Y], coordinates of the points in 2D
% Z: z data set.

a = x(1);
b = x(2);
c = x(3);
d = x(4);
x0 = x(5);
y0 = x(6);

X = XY(:,1);
Y = XY(:,2);

val = exp(-(a*(X-x0).^2 + 2*b*(X-x0).*(Y-y0) + c*(Y-y0).^2)) + d;

end % fitGauss2d