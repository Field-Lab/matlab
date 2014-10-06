function r=ator(a,b)

%ATOR	Inverse Levinson-Durbin recursion.

%----

%USAGE: r=ator(a,b)

%

%	Finds the autocorrelation sequence r(k) of an

%       autoregressive process that is generated by 

%	filtering unit variance white noise with the filter

%		H(z)=b(0)/A(z)

%       With r = ator(a) the autocorrelation sequence is 

%	normalized so that r(0)=1.

%

%  see also ATOG, GTOA, GTOR, RTOA, RTOG

%

%---------------------------------------------------------------

% copyright 1996, by M.H. Hayes.  For use with the book 

% "Statistical Digital Signal Processing and Modeling"

% (John Wiley & Sons, 1996).

%---------------------------------------------------------------



   p=length(a)-1;

   gamma=atog(a);

   r=gtor(gamma);

   if nargin == 2,

   r = r*b^2/prod(1-abs(gamma).^2);

   end;

   

