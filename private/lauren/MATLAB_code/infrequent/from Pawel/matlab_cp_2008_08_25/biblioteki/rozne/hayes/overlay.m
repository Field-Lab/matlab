function Px=overlay(N,omega,A,sigma,num)

%OVERLAY Forms multiple power spectral density estimates.  

%-------

%USAGE	Px=overlay(N,omega,A,sigma,num)

%

%	Calculates the periodogram using an ensemble of realizations 

%	of a process consisting of a sum of complex exponential in 

%	white noise.

%

%          N     :  length of the signal

%          omega :  vector containing the sinusoid frequencies

%          A     :  vector containing the sinusoid frequencies

%          sigma :  the variance of the white noise

%          num   :  size of the ensemble (number of periodograms)

%

%---------------------------------------------------------------

% copyright 1996, by M.H. Hayes.  For use with the book 

% "Statistical Digital Signal Processing and Modeling"

% (John Wiley & Sons, 1996).

%---------------------------------------------------------------



jj=length(omega);

n=1:N;

for i=1:num

     x=sigma*randn(1,N);

     for j=1:jj

           phi=2*pi*rand(1);

           x=x+A(j)*sin(omega(j)*n+phi);

           end;

     Px(:,i)=per(x);   % This command may be replaced 

                       % with another m-file.

     end;

