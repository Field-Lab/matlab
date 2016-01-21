% Function that does ideal sinc interpolation
% Chaitu Ekanadham

% Arguments:
% f - the interpolator filter with length m
% timepts - n points at which samples are taken (indices between 1 and m)

% Returns:
% A - the m x n interpolation matrix. i.e. if x is a column vector of sample
% values at timepts, then A*x is the (bandlimited) function at taxis,
% assuming zero outside.

function A = interpmtx(f,timepts,T)


% f should be odd length
m = length(f);
n = length(timepts);
flong = [f(:); zeros(T,1)]; % pad with trailing zeros
midrange = floor(m/2):floor(m/2)+T-1;
1;

% For i=1,...,n, the ith col of A will be f(k-i):1 <= k <= m assuming zero padding
A = zeros(T,n);
for i=1:n
    x = circshift(flong,timepts(i)-1);
    A(:,i) = x(midrange);
end