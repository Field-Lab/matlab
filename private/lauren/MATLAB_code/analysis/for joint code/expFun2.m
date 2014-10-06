function error = expFun2(nlParams, xValues, yValues)

% arguments:
%   a: exponential coefficient to fit
%
%   data: 2-d array, with data(1,:) vector of x-values and data(2,:) vector of y-values
%
% returns: distance between vector of data values and vector of calculated values based on current
% params

% linear regression to find linear parameters
X = [ones(size(xValues))' exp(nlParams(1)*xValues)' exp(nlParams(2)*xValues)'];

linParams = X\yValues;

a = nlParams(1);
b = nlParams(2);

c = linParams(1);
d = linParams(2);
e = linParams(3);


% plug in b
p = c + d*exp(a*xValues) + e*exp(b*xValues);

% clf
% hold on
% plot(xValues, p, 'k')
% plot(xValues, yValues,'b')
% hold off

error = norm(p - yValues');