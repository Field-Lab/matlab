function [a,b,err] = shanks(x,p,q)

%SHANKS	Model a signal using Shanks' method 

%----

%Usage: [a,b,err] = shanks(x,p,q) 

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

%  see also ACM, COVM, PADE, PRONY

%

%---------------------------------------------------------------

% copyright 1996, by M.H. Hayes.  For use with the book 

% "Statistical Digital Signal Processing and Modeling"

% (John Wiley & Sons, 1996).

%---------------------------------------------------------------



x   = x(:);

N   = length(x);

if p+q>=length(x), error('Model order too large'), end

a   = prony(x,p,q);

u   = [1; zeros(N-1,1)];

g   = filter(1,a,u);

G   = convm(g,q+1);

b   = G(1:N,:)\x;

err = x'*x-x'*G(1:N,1:q+1)*b;

end;



