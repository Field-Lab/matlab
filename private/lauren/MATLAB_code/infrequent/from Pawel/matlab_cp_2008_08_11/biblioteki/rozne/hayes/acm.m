function [a,err] = acm(x,p)

%ACM	Find an all-pole model using the autocorrelation method

%----

%Usage: [a,err] = acm(x,p) 

%

%	The input sequence x is modeled as the unit sample response of

%	a filter having a system function of the form

%		H(z) = b(0)/A(z) 

%	where the coefficients of A(z) are contained in the vector

%		a=[1, a(1), ... a(p)]

%	The input p defines the number of poles in the model.

%	The modeling error is returned in err.

%	The numerator b(0) is typically set equal to the square

%	root of err.

%

%  see also COVM, PADE, PRONY, SHANKS

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

Xq  = X(1:N+p-1,1:p);

a   = [1;-Xq\X(2:N+p,1)];

err = abs(X(1:N+p,1)'*X*a);

end;



