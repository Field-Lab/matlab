function [gamma,err] = fcov(x,p)

%FCOV	Forward covariance algorithm.

%----

%USAGE	[gamma,err] = fcov(x,p) 

%

%	An all-pole model of order p for the input sequence x is

%	found using the forward covariance method.  The outputs

%	are the reflection coefficients

%		gamma = [g(1), g(2), ... , g(p)]

%	corresponding to the all-pole model, and the model error

%	which is returned in err.

%

%  see also COVM, MCOV and BURG

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

N=N-1;

for j=1:p;

      gamma(j) = -eminus'*eplus/(eminus'*eminus);

      temp1    =  eplus  + gamma(j)*eminus;

      temp2    =  eminus + conj(gamma(j))*eplus;

      err(j)   =  temp1'*temp1;

      eplus    =  temp1(2:N);

      eminus   =  temp2(1:N-1);

      N=N-1;

      end;

