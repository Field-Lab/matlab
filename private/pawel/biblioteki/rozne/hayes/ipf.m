function [a,b,err] = ipf(x,p,q,n,a)

%IPF	Pole-zero signal modeling using iterative prefiltering.

%---

%USAGE	[a,b,err] = ipf(x,p,q,n,a) 

%

%	The sequence x is modeled as the unit sample response of

%	a filter having a system function of the form

%		H(z) = B(z)/A(z) 

%	The polynomials B(z) and A(z) are formed from the vectors

%		b=[b(0), b(1), ... b(q)]

%		a=[1   , a(1), ... a(p)]

%	The inputs are

%          x   :  input sequence to be modeled

%          p   :  number of poles in the model

%          q   :  number of zeros in the model

%          n   :  number of iterations to be used to find A(z) and B(z)	

%          a   :  initial estimate for the denominator (default is zero).

%	and the outputs are

%          a   :  vector of coefficients for A(z)

%          b   :  vector of coefficients for B(z)

%          err :  model error

%

%---------------------------------------------------------------

% copyright 1996, by M.H. Hayes.  For use with the book 

% "Statistical Digital Signal Processing and Modeling"

% (John Wiley & Sons, 1996).

%---------------------------------------------------------------



x   = x(:);

N   = length(x);

if p+q>=length(x), error('Model order too large'), end

if nargin < 5 

    a   = prony(x,p,q);

    end

delta = [1; zeros(N-1,1)];



for i=1:n

    f   = filter(1,a,x);

    g   = filter(1,a,delta);

    u   = convm(f,p+1); 

    v   = convm(g,q+1);

    ab  = -[u(1:N,2:p+1) -v(1:N,:) ]\u(1:N,1);

    a   = [1; ab(1:p)];

    b   = ab(p+1:p+q+1);

    err = norm( u(1:N,1) + [u(1:N,2:p+1) -v(1:N,:)]*ab);

    end;







