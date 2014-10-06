function Px = bt_pc(x,p,M)

%BT_PC	Frequency estimation using principal components Blackman-Tukey.

%-----

%USAGE	Px = bt_pc(x,p,M)

%

%	The spectrum of a process x is estimated using a principal

%	components analysis of the autocorrelation matrix.

%	The model for the process is that x(n) consists of a sum of 

%	complex exponentials in white noise.  

%	After a principle components analysis, the principal eigenvectors

%	are used in the Blackman-Tukey estimate.

%

%          x  : input sequence

%          p  : number of complex exponentials in x

%          M  : size of autocorrelation matrix

%

%	The spectrum estimate is returned in Px using a dB scale.

%

%---------------------------------------------------------------

% copyright 1996, by M.H. Hayes.  For use with the book 

% "Statistical Digital Signal Processing and Modeling"

% (John Wiley & Sons, 1996).

%---------------------------------------------------------------



   x   = x(:);

   if M<p+1, error('Specified size of R is too small'), end

   R=covar(x,M);

   [v,d]=eig(R);

   [y,i]=sort(diag(d));

   Px=0;

   for j=M-p+1,M;

       Px=Px+abs(fft(v(:,i(j)),1024))*sqrt(real(y(j)));

       end;

   Px=20*log10(Px)-10*log10(M);

   end;

































