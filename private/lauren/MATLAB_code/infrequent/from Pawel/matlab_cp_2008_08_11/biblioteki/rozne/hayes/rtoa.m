function [a,epsilon]=rtoa(r)

%RTOA	Levinson-Durbin recursion.

%----

%USAGE:	[a,epsilon]=rtoa(r)

%

%	Solves the Toeplitz normal equations

%		R a = epsilon [1 0 ... 0]'

% 	where R=toeplitz(r) is a Toeplitz matrix that contains

%	the autocorrelation sequence r(k).  The ouputs are the

%	coefficients of the all-pole model a(k) and the constant

%	epsilon.

%

%  see also ATOG, ATOR, GTOA, GTOR, RTOG

%

%---------------------------------------------------------------

% copyright 1996, by M.H. Hayes.  For use with the book 

% "Statistical Digital Signal Processing and Modeling"

% (John Wiley & Sons, 1996).

%---------------------------------------------------------------



r=r(:);

p=length(r)-1;

a=1;

epsilon=r(1);

for j=2:p+1;

	gamma=-r(2:j)'*flipud(a)/epsilon;

	a=[a;0] + gamma*[0;conj(flipud(a))];

	epsilon=epsilon*(1 - abs(gamma)^2);

	end











