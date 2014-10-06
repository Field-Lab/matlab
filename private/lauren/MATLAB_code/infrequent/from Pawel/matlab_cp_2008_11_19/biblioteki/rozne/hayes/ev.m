function Px = ev(x,p,M)

%EV	Frequency estimation using the eigenvector method

%--

%USAGE	Px = ev(x,p,M)

%

%	The input sequence x is assumed to consist of p complex

%	exponentials in white noise.  The frequencies of the

%	complex exponentials and the variance of the white noise

%	are estimated using the eigenvector method.  

%

%	x : input sequence

%	p : Number of complex exponential in x

%	M : Size of the autocorrelation matrix to use in

%	      estimating the complex exponential frequencies

%

%	The frequency estimates are found from the peaks of the

%	pseudospectrum Px.

%

%  see also PHD, MUSIC, and MIN_NORM

%

%---------------------------------------------------------------

% copyright 1996, by M.H. Hayes.  For use with the book 

% "Statistical Digital Signal Processing and Modeling"

% (John Wiley & Sons, 1996).

%---------------------------------------------------------------



   x   = x(:);

   if N<p+1, error('Specified size of R is too small'), end

   R = covar(x,N);

   [v,d]=eig(R);

   [y,i]=sort(diag(d));

   Px=0;

   for j=1:N-p

       Px=Px+abs(fft(v(:,i(j)),1024)).^2/abs(y(j));

       end;

   Px=-10*log10(Px);

   end;















