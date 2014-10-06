function [error projection] = expFitter2(trace)

nlParamsStart(1) = -0.5; %initial guess for nonlinear decay parameter
nlParamsStart(2) = -0.1; %initial guess for nonlinear increase parameter

xValues = 1:length(trace);
yValues = trace;

warning off MATLAB:rankDeficientMatrix %to suppress warning while searching for best nonlinear parameter
estimatedNlParams = fminsearch(@(nlParams)expFun2(nlParams, xValues, yValues), nlParamsStart);
warning on MATLAB:rankDeficientMatrix

a = estimatedNlParams(1);
b = estimatedNlParams(2);

X = [ones(size(xValues))' exp(a*xValues)' exp(b*xValues)'];

lParams = X\yValues;

projection = lParams(1) + lParams(2)*exp(estimatedNlParams(1)*xValues) + lParams(3)*exp(estimatedNlParams(2)*xValues);
error = norm(projection-trace');
% 
% figure(2)
% hold on
% plot(xValues, projection, 'k')
% plot(xValues, yValues, 'b')
% hold off
% 
% error2 = expFun2(estimatedNlParams, xValues, yValues);
% disp()

%disp(['a=' num2str(estimatedNlParams(1)) ', b=' num2str(estimatedNlParams(2))])