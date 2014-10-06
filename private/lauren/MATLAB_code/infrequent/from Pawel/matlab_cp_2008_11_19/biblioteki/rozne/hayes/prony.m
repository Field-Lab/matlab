function [a,b,err] = prony(x,p,q)

%PRONY	Model a signal using Prony's method

%----

%Usage: [a,b,err] = prony(x,p,q) 

%

%	The sequence x is modeled as the unit sample response of

%	a filter having a system function of the form

%		H(z) = B(z)/A(z) 

%	The polynomials B(z) and A(z) are formed from the vectors

%		b=[b(0), b(1), ... b(q)]

%		a=[1   , a(1), ... a(p)]

%	The input q defines the number of zeros in the model

%	and p defines the number of poles.  The modeling error

%	is returned in err.

%

%  see also ACM, COVM, PADE, SHANKS

%

%---------------------------------------------------------------

% copyright 1996, by M.H. Hayes.  For use with the book 

% "Statistical Digital Signal Processing and Modeling"

% (John Wiley & Sons, 1996).

%---------------------------------------------------------------



x   = x(:);

N   = length(x);

if p+q>=length(x), error('Model order too large'), end

X   = convm(x,p+1);

Xq  = X(q+1:N+p-1,1:p);

a   = [1;-Xq\X(q+2:N+p,1)];

b   = X(1:q+1,1:p+1)*a;

err = x(q+2:N)'*X(q+2:N,1:p+1)*a;

end;









