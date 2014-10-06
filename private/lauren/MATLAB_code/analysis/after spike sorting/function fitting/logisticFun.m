function error = logisticFun(params, data)

% arguments:
%   params: parameters a, b, and c to be fit
%
%   data: 2-d array, with data(1,:) vector of stimulus amplitudes and data(2,:) vector of success
%   rates
%
% returns: distance between vector of data values and vector of calculated values based on current
% params


lam = params(1);
th = params(2);



p = zeros(1, size(data,2));

for i = 1:size(data,2)
    if (data(1,i)-th)/lam >0
    %p(i) = 1./(1 + exp(-(((data(1,i)-th)./lam).^2)));
        p(i) = 1./(1 + exp(-(((data(1,i)-th)./lam).^2)));
    else
        p(i) = 1./(1 + exp((((data(1,i)-th)./lam).^2)));
    end
end

%figure
%plot(data(1,:), p)

error = norm(p - data(2,:));