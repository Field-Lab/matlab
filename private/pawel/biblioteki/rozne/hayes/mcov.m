function [a,err] = mcov(x,p)

%MCOV	All-pole signal modeling using the modified covariance method.

%----

%USAGE	[a,err] = mcov(x,p) 

%

%	An all-pole of order p is found for the input sequence x using

%	the modified covariance method.  The model is of the form

%		H(z) = b(0)/A(z) 

%	The coefficients of A(z) are returned in the vector

%		a=[1, a(1), ... a(p)]

%	and the modeling error is returned in err.

%

%  See also:  FCOV, COVM, and BURG

%

%---------------------------------------------------------------

% copyright 1996, by M.H. Hayes.  For use with the book 

% "Statistical Digital Signal Processing and Modeling"

% (John Wiley & Sons, 1996).

%---------------------------------------------------------------

%



x   = x(:);

N   = length(x);

if p>=length(x), error('Model order too large'), end

X   = toeplitz(x(p+1:N),flipud(x(1:p+1)));

R   = X'*X;

R1  = R(2:p+1,2:p+1);

R2  = flipud(fliplr(R(1:p,1:p)));

b1  = R(2:p+1,1);

b2  = flipud(R(1:p,p+1));

a   = [1 ; -(R1+R2)\(b1+b2)];

err = R(1,:)*a+fliplr(R(p+1,:))*a;

end;







