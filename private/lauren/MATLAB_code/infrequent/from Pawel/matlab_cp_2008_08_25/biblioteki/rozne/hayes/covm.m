function [a,err] = covm(x,p)

%COVM	Find an all-pole model using the covariance method

%----

%USAGE:	[a,err] = covm(x,p) 

%

%	An all-pole of order p is found for the input sequence 

%	x using the covariance method.  The model is of the form

%		H(z) = b(0)/A(z) 

%	The coefficients of A(z) are returned in the vector

%		a=[1, a(1), ... a(p)]

%	and the modeling error is returned in err.

%

%  see also FCOV, MCOV, and BURG

%

%---------------------------------------------------------------

% copyright 1996, by M.H. Hayes.  For use with the book 

% "Statistical Digital Signal Processing and Modeling"

% (John Wiley & Sons, 1996).

%---------------------------------------------------------------



x   = x(:);

N   = length(x);

if p>=length(x), error('Model order too large'), end

X   = convm(x,p+1);

Xq  = X(p:N-1,1:p);

a   = [1;-Xq\X(p+1:N,1)];

err = abs(X(p+1:N,1)'*X(p+1:N,:)*a);

end;

