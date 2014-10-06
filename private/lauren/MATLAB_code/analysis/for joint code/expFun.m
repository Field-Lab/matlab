function error = expFun(a, xValues, yValues)

% arguments:
%   a: exponential coefficient to fit
%
%   data: 2-d array, with data(1,:) vector of x-values and data(2,:) vector of y-values
%
% returns: distance between vector of data values and vector of calculated values based on current
% params

% linear regression to find linear parameters
X = [ones(size(xValues))' exp(a*xValues)'];

linParams = X\yValues;
    
b = linParams(1);
c = linParams(2);

% plug in b
p = b + c*exp(a*xValues);

error = norm(p - yValues');