function [slope y_0] = fitLine(xVals, yVals)

%% linear regression of x and y values to give slope and y-intercept (y_0)

A = ones(length(xVals), 2);
A(:,1) = xVals;
if size(yVals, 2) == 1
    b = yVals;
else
    b = yVals';
end

x = A\b;

slope = x(1);
y_0 = x(2);
