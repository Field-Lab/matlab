function [wspol,r] = filtr_faza_f(b,ilosc);

a=size(b);
if a(1)>1 
   b=b';
end

opoznienie=(a(2)-1)/2;

x=zeros(1,ilosc);
x(1,1:a(2))=b;

w=fft(x);

k=[0:(ilosc-1)];
%w=w.*exp(2*pi*j*opoznienie/ilosc*k);

wspol=w;
r=ifft(w);