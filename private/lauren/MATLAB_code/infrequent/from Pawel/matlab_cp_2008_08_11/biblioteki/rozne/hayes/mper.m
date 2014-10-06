function Px = mper(x,win,n1,n2)

%MPER	Spectrum estimation using the modified periodogram.

%----

%USAGE	Px = mper(x,win,n1,n2) 

%

%	The spectrum of a process x is estimated using the modified

%	periodogrm.

%

%	x   : input sequence

%       n1  : starting index, x(n1)

%	n2  : ending index, x(n2)

%	win : The window type 

%       	1 = Rectangular

%		2 = Hamming

%		3 = Hanning

%		4 = Bartlett

%		5 = Blackman

%

%       If n1 and n2 are not specified the periodogram of the entire

%       sequence is computed.

%

%	The modified periodogram is returned in Px using a linear scale.

%

%  see also BART, PER, SPER, and WELCH

%

%---------------------------------------------------------------

% copyright 1996, by M.H. Hayes.  For use with the book 

% "Statistical Digital Signal Processing and Modeling"

% (John Wiley & Sons, 1996).

%---------------------------------------------------------------



x   = x(:);

if nargin == 2

    n1 = 1;  n2 = length(x);  end;

N  = n2 - n1 +1;

w  = ones(N,1);

if (win == 2) w = hamming(N);

   elseif (win == 3) w = hanning(N);

   elseif (win == 4) w = bartlett(N);

   elseif (win == 5) w = blackman(N); 

   end;

U  = norm(w)^2/N;

xw = x(n1:n2).*w;

Px = abs(fft(xw,1024)).^2/(N*U);

Px(1)=Px(2);

end;



