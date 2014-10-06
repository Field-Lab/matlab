function Px = per(x,n1,n2)

%PER	Estimate the spectrum of a process using the periodogram

%-----------

%USAGE	Px = per(x,n1,n2) 

%

%	The spectrum of a process x is estimated using the periodogram.

%	

%	x   :  input sequence

%       n1  :  starting index, x(n1)

%	n2  :  ending index, x(n2)

%

%	The periodogram is returned in Px using a linear scale.

%

%       If n1 and n2 are not specified the periodogram of the entire

%       sequence is computed.

%

%  see also BART, MPER, WELCH, and SPER

%

%---------------------------------------------------------------

% copyright 1996, by M.H. Hayes.  For use with the book 

% "Statistical Digital Signal Processing and Modeling"

% (John Wiley & Sons, 1996).

%---------------------------------------------------------------



x   = x(:);

if nargin == 1

    n1 = 1;  n2 = length(x);  end;

    Px = abs(fft(x(n1:n2),1024)).^2/(n2-n1+1);

    Px(1)=Px(2);

end;











