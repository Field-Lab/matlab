function Px = music(x,p,M)

%MUSIC	Frequency estimation using the MUSIC algorithm.

%-----

%USAGE	Px=music(x,p,M)

%

%	The input sequence x is assumed to consist of p complex

%	exponentials in white noise.  The frequencies of the

%	complex exponentials and the variance of the white noise

%	are estimated using the MUSIC algorithm.  

%

%	x : input sequence

%       p : number of complex exponentials to find

%	M : number of noise eigenvectors to use

%

%	The frequency estimates are found from the peaks of the

%	pseudospectrum Px.

%

%  see also PHD, EV, and MIN_NORM

%

%---------------------------------------------------------------

% copyright 1996, by M.H. Hayes.  For use with the book 

% "Statistical Digital Signal Processing and Modeling"

% (John Wiley & Sons, 1996).

%---------------------------------------------------------------



   x   = x(:);

   if M<p+1 | length(x)<M, error('Size of R is inappropriate'), end

   R = covar(x,M);

   [v,d]=eig(R);

   [y,i]=sort(diag(d));

   Px=0;

   for j=1:M-p

       Px=Px+abs(fft(v(:,i(j)),1024));

       end;

   Px=-20*log10(Px);

   end;









