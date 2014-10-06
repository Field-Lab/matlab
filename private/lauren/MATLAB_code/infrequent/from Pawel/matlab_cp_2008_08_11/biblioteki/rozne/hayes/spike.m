function [h,err] = spike(g,n0,n)

%SPIKE	Finds the FIR least squares inverse of g(n)

%----

%USAGE:	[h,err] = spike(g,n0,n) 

%

%	g   : the filter that is to be inverted

%	n0  : the delay used in the minimization of

%		e(n)=g(n)*g(n) - delta(n-n0)

%	n   : order (number of coefficients) in h

%	h   : the FIR least squares inverse of g

%	err : squared error in the approximation

%

%---------------------------------------------------------------

% copyright 1996, by M.H. Hayes.  For use with the book 

% "Statistical Digital Signal Processing and Modeling"

% (John Wiley & Sons, 1996).

%---------------------------------------------------------------



g   = g(:);

m   = length(g);

if m+n-1<=n0, error('Delay too large'), end

G   = convm(g,n);

d   = zeros(m+n-1,1);

d(n0+1) = 1;

h   = G\d;

err = 1 - G(n0+1,:)*h;

end;

