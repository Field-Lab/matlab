function [gamma,epsilon] = rtog(r)

%RTOG	Levinson-Durbin recursion.

%----

%USAGE:	[gamma,epsilon] = rtog(r)

%

%	Solves the Toeplitz normal equations

%		R a = epsilon [1 0 ... 0]'

% 	where R=toeplitz(r) is a Toeplitz matrix that contains

%	the autocorrelation sequence r(k), and then maps the

%	filter coefficients into the corresponding reflection

%	coefficient sequence gamma.

%

%  see also ATOG, ATOR, GTOA, GTOR, RTOA

%

%---------------------------------------------------------------

% copyright 1996, by M.H. Hayes.  For use with the book 

% "Statistical Digital Signal Processing and Modeling"

% (John Wiley & Sons, 1996).

%---------------------------------------------------------------



   [a,epsilon]=rtoa(r);

   gamma=atog(a);

   end; 





