function [gamma,err] = burg(x,p)

%BURG	All-pole modeling using the Burg algorithm.

%----

%USAGE	[gamma,err] = burg(x,p) 

%

%	An all-pole of order p is found for the input sequence 

%	x using the Burg algorithm.  The outputs are the reflection

%	coefficients

%		gamma = [g(1), g(2), ... , g(p)]

%	corresponding to an all-pole model, and the Burg error, 

%	which is returned in err.

%

%  see also COVM, FCOV and MCOV 

%

%---------------------------------------------------------------

% copyright 1996, by M.H. Hayes.  For use with the book 

% "Statistical Digital Signal Processing and Modeling"

% (John Wiley & Sons, 1996).

%---------------------------------------------------------------



x = x(:);

N=length(x);

eplus  = x(2:N);

eminus = x(1:N-1);

N=length(eplus);

for j=1:p;

      gamma(j) = -2*eplus'*eminus/(eplus'*eplus+eminus'*eminus);

      temp1    = eplus  + gamma(j)*eminus;

      temp2    = eminus + gamma(j)*eplus;

      error(j) = temp1'*temp1+temp2'*temp2;

      eplus    = temp1(2:N);

      eminus   = temp2(1:N-1);

      N=N-1;

      end;

