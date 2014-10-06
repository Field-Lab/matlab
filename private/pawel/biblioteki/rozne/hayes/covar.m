function R=covar(x,p)

%COVAR	Generates a covariance matrix/

%-----

%USAGE	R=covar(x,p)

%

%	Generates a p x p covariance matrix for the sequence x.

%

%---------------------------------------------------------------

% copyright 1996, by M.H. Hayes.  For use with the book 

% "Statistical Digital Signal Processing and Modeling"

% (John Wiley & Sons, 1996).

%---------------------------------------------------------------



x = x(:);

   m = length(x);

   x = x - ones(m,1)*(sum(x)/m);

   R = convm(x,p)'*convm(x,p)/(m-1);

   end;

