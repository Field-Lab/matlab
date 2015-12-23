function error = weibullFun(params, data, fixedParam, fixedParamID)

% arguments:
%   params: parameters a, b, and c to be fit
%
%   data: 2-d array, with data(1,:) vector of stimulus amplitudes and data(2,:) vector of success
%   rates
%
% returns: distance between vector of data values and vector of calculated values based on current
% params

if nargin > 2
    if fixedParamID == 1
        k = fixedParam;
        lam = params(1);
        th = params(2);
    elseif fixedParamID == 2
        k = params(1);
        lam = fixedParam;
        th = params(2);
    else
        k = params(1);
        lam = params(2);
        th = fixedParam;
    end
else
    k = params(1);
    lam = params(2);
    th = params(3);
end

p = zeros(1, size(data,2));

for i = 1:size(data,2)
    if (data(1,i) - th)/lam >= 0
        p(i) = 1 - exp(-((data(1,i)-th)/lam)^k);
    else
        p(i) = 0;
    end
end

%figure
%plot(data(1,:), p)

error = norm(p - data(2,:));