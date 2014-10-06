function val = boxsplinen(x, y, n)
% BOXSPLINEN
% usage: val = boxsplinen(x, y, n)
%
% X and Y are equal size coordinate vectors/matrices, N is the spline order.  
% For N == 2, there is a much faster C++ implementation available in 
% BOXSPLINE2 and BOXSPLINE2_PARALLEL.
%
% Slightly cleaned and optimized copy from:
%   Three-Directional Box-Splines: Characterization and Efficient
%   Evaluation, Laurent Condat and Dimitri Van De Ville, IEEE Signal
%   Processing Letters, 2006
%
% See also: BOXSPLINE2, BOXSPLINE2_PARALLEL
%
% 2011-08, phli
%

x = -abs(x);
y =  abs(y);

u = x - y/sqrt(3);
v = x + y/sqrt(3);

id = v > 0;
v(id) = -v(id);
u(id) = u(id) + v(id);

id = v > u/2; 
v(id) = u(id) - v(id);


% Only calculate for the area within the support
within = v <= n & v >= -n;
U = u(within);
V = v(within);

val = zeros(size(x));
val(within) = boxsplinen_(U, V, n);



function val = boxsplinen_(u, v, n)

% Precalculate and store in memory...
fact = factorial(0:(3*n-2)); % shifted by 1, i.e. fact(1) = factorial(0)

val = zeros(size(u));
for K = -n:(ceil(max(max(u)))-1)

    for L = -n:(ceil(max(max(v)))-1)

        aux  = abs(v-L-u+K);
        aux2 = (u-K+v-L-aux)/2;
        aux2(aux2<0) = 0;
        
        for i = 0:min(n+K,n+L)
            
            coeff = (-1)^(K+L+i) * nchoosek(n,i-K) * nchoosek(n,i-L) * nchoosek(n,i);            
            for d = 0:n-1
                val = val + coeff*nchoosek(n-1+d,d)/ ...
                    fact(2*n+d)/fact(n-d)* ...
                    aux.^(n-1-d).*aux2.^(2*n-1+d);
            end
            
        end
        
    end
end