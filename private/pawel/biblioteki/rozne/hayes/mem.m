function Px = mem(x,p)

%MEM	Spectrum estimation using the Maximum Entropy Method (MEM).

%---

%USAGE	Px = mem(x,p) 

%

%	The spectrum of a process x is estimated using the maximum

%	entropy method, which uses the autocorrelation method to 

%	find a pth-order all-pole model for x(n), and then forms 

%	the estimate of the spectrum as follows:

%		Px = b^2(0)/|A(omega)|^2 

%	The spectrum estimate is returned in Px using a dB scale.

%          x  :  Input sequence

%          p  :  Order of the all-pole model

%

%---------------------------------------------------------------

% copyright 1996, by M.H. Hayes.  For use with the book 

% "Statistical Digital Signal Processing and Modeling"

% (John Wiley & Sons, 1996).

%---------------------------------------------------------------



[a,e] = acm(x,p);

   Px = 10*log10(e/length(x))-20*log10(abs(fft(a,1024)));

end;









