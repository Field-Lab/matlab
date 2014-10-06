function [y] = filtrhigh2(s,t,f0);
%FILTRHIGH2 Summary of this function goes here
%   s - sygna³
%   t - wektor czasu
%   f0 - czêstotliwoœæ graniczna
F=0:1/10000:1;
H2=[zeros(1,f0+1),ones(1,10000-f0)];
h2=fir2(200,F,H2);
%figure
%freqz(h2,1)
tic
y=fftfilt(h2,s);
toc
figure
plot (t,y)
end

