function error = expFitter(trace)

aGuess = -1; %initial guess for nonlinear parameter
params = zeros(1,3);

xValues = 1:length(trace);
yValues = trace;

warning off MATLAB:rankDeficientMatrix %to suppress warning while searching for best nonlinear parameter
params(1) = fminsearch(@(nlParam)expFun(nlParam, xValues, yValues), aGuess);
warning on MATLAB:rankDeficientMatrix

X = [ones(size(xValues))' exp(params(1)*xValues)'];
params(2:3) = X\yValues;
projection = params(2) + params(3)*exp(params(1)*xValues);
error = norm(projection-trace');