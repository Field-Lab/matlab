function Px = welch(x,L,over,win)

%WELCH	Spectrum estimation using Welch's method.

%----

%USAGE	Px = welch(x,L,over,win)

%

%	The spectrum of a process x is estimated using Welch's method

%	of averaging modified periodograms.

%

%	x   : input sequence

%	L   : section length 

%	over: amount of overlap, where 0<over<1, 

%	win : The window type 

%       	1 = Rectangular

%		2 = Hamming

%		3 = Hanning

%		4 = Bartlett

%		5 = Blackman

%

%	Welch's estimate of the power spectrum is returned in Px 

%	using a linear scale.

%

%  see also BART, MPER, PER, and SPER

%

%---------------------------------------------------------------

% copyright 1996, by M.H. Hayes.  For use with the book 

% "Statistical Digital Signal Processing and Modeling"

% (John Wiley & Sons, 1996).

%---------------------------------------------------------------



if (nargin <= 3) win=1; end;

if (nargin <= 2) over=0; end;

if (nargin == 1) L=length(x); end

if (over >= 1) | (over < 0)  

   error('Overlap is invalid'), end

n1 = 1;

n2 = L;

n0 = (1-over)*L;

nsect=1+floor((length(x)-L)/(n0));

Px=0;

for i=1:nsect

    Px = Px + mper(x,win,n1,n2)/nsect;

    n1 = n1 + n0;  

    n2 = n2 + n0;

    end;









