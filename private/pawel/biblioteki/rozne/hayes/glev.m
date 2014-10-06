function x=glev(r,b)

%GLEV  	Levinson recursion.

%----

%USAGE: x = glev(r,b)

%

%	The Levinson recursion solves the Toeplitz equations

%		R x = b

%	where R=toeplitz(r) and b is an arbitrary vector.

%	

%  see also ATOG, ATOR, GTOA, GTOR, RTOA, RTOG

%

%---------------------------------------------------------------

% copyright 1996, by M.H. Hayes.  For use with the book 

% "Statistical Digital Signal Processing and Modeling"

% (John Wiley & Sons, 1996).

%---------------------------------------------------------------



r=r(:);

p=length(b);

a=1;

x=b(1)/r(1);

epsilon=r(1);

for j=2:p;

	g=r(2:j)'*flipud(a)

	gamma(j-1)=-g/epsilon

	a=[a;0] + gamma(j-1)*[0;conj(flipud(a))]

	epsilon=epsilon*(1 - abs(gamma(j-1))^2)

	delta=r(2:j)'*flipud(x)

	q=(b(j)-delta)/epsilon

	x=[x;0] + q*[conj(flipud(a))]

	end

