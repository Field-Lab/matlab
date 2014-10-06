function a=gtoa(gamma)

%GTOA	Step-up recursion

%----

%USAGE	a = gtoa(gamma) 

%

%	The direct-form filter coefficients

%		a=[1 a(1) ... a(p)]

%	are derived from the reflection coefficients gamma

%	using the step-up recursion.

%

%  see also ATOG, ATOR, GTOR, RTOA, RTOG

%

%---------------------------------------------------------------

% copyright 1996, by M.H. Hayes.  For use with the book 

% "Statistical Digital Signal Processing and Modeling"

% (John Wiley & Sons, 1996).

%---------------------------------------------------------------



a=1;

gamma=gamma(:);

p=length(gamma);

for j=2:p+1;

	a=[a;0] + gamma(j-1)*[0;conj(flipud(a))];

	end



