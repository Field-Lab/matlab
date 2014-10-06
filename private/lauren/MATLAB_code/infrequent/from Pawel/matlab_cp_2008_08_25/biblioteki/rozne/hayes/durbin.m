function b = durbin(x,p,q)

%DURBIN	Find a moving average model using Durbin's method

%----

%USAGE:	b = durbin(x,p,q) 

%

%	A moving average model for the input sequence x is found

%	using Durbin's method.  The coefficients are returned in

%	the vector

%		b = [b(0), b(1), ... , b(q) ]

%	

%	q :  order of the moving average model

%	p :  order of the all-pole model used to approximation 1/B(z).

%	     This parameter should be at least 4*q.

%

%---------------------------------------------------------------

% copyright 1996, by M.H. Hayes.  For use with the book 

% "Statistical Digital Signal Processing and Modeling"

% (John Wiley & Sons, 1996).

%---------------------------------------------------------------



x   = x(:);

if p>=length(x), error('Model order too large'), end

[a,epsilon] = acm(x,p);

[b,epsilon] = acm(length(x)*a/sqrt(epsilon),q);

b = b*length(x)/sqrt(epsilon);

end;



