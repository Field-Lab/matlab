function Px = min_norm(x,p,M)

%MIN_NORM Frequency estimation using the minimum norm algorithm.

%--------

%USAGE	Px = min_norm(x,p,M)

%

%	The input sequence x is assumed to consist of p complex

%	exponentials in white noise.  The frequencies of the

%	complex exponentials and the variance of the white noise

%	are estimated using the minimum norm algorithm.  

%

%	x : input sequence

%	p : Number of complex exponential in x

%	M : Size of the autocorrelation matrix to use in

%	      estimating the complex exponential frequencies

%

%	The frequency estimates are found from the peaks of the

%	pseudospectrum Px.

%

%  see also PHD, EV, and MUSIC

%

%---------------------------------------------------------------

% copyright 1996, by M.H. Hayes.  For use with the book 

% "Statistical Digital Signal Processing and Modeling"

% (John Wiley & Sons, 1996).

%---------------------------------------------------------------

%



   x   = x(:);

   if N<p+1, error('Specified size of R is too small'), end

   R=covar(x,N);

   [v,d]=eig(R);

   [y,i]=sort(diag(d));

   for j=1:N-p

       V=[V,v(:,i(j))];

       end;

   a=V*V(1,:)';

   Px=-20*log10(abs(fft(a,1024)));

   end;

