function error = erfFunPlotIter(params, data)

% arguments:
%   params: parameters a and b to be fit
%
%   data: 2-d array, with data(1,:) vector of stimulus amplitudes and data(2,:) vector of success
%   rates
%
% returns: distance between vector of data values and vector of calculated values based on current
% params

a = params(1);
b = params(2);

p = zeros(1, size(data,2));
for i = 1:size(data,2)
    p(i) = 0.5 + 0.5*erf(a*data(1,i) + b);
end

error = norm(p - data(2,:));

clf
hold on
plot(data(1,:), data(2,:), 'k.')

xProj = data(1,1):0.001:data(1,end);
pFull = 0.5*(1 + erf(a*xProj + b));

plot(xProj, pFull, 'b')

hold off
drawnow