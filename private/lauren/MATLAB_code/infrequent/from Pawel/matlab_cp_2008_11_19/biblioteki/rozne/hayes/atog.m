function gamma=atog(a)

%ATOG	Step-down recursion

%----

%USAGE: gamma=atog(a)

%

%	Finds the reflection coefficients gamma from the

%	direct-form filter coefficients a(k).

%

%  see also ATOR, GTOA, GTOR, RTOA, RTOG

%

%---------------------------------------------------------------

% copyright 1996, by M.H. Hayes.  For use with the book 

% "Statistical Digital Signal Processing and Modeling"

% (John Wiley & Sons, 1996).

%---------------------------------------------------------------



a=a(:);

p=length(a);

a=a(2:p)/a(1);

gamma(p-1)=a(p-1);

for j=p-1:-1:2;

	a=(a(1:j-1) - gamma(j)*flipud(conj(a(1:j-1))))./ ...

	  (1 - abs(gamma(j))^2);

	gamma(j-1)=a(j-1);

	end

