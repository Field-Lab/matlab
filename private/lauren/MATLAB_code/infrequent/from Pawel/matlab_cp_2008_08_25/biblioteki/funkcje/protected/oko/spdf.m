function [f,w]=spdf(x,N,M,fs,lsb);
%Funkcja liczaca gestosc widmowa szumow [V*V/Hz];
%x - probki;
%N - ilosc czestotliwosci;
%M - krok, o jaki przesuwa sie okno czasowe;
%fp - czestotliwosc probkowania;
%lsb - stala normalizacyjna wyrazona w voltach
%wyjscie:
%f - czestotliwosci
%w - wartosci gestosci widmowej
N;
M;
fs;

a=size(x);
if a(1)>1 
   x=x';
end
a=size(x);
length=a(2);

steps=floor((length-N)/M)+1;
df=fs/N;
f=[0:(N-1)]/N*fs;
w=zeros(1,N);

for i=1:steps
   i1=1+(i-1)*M;
   i2=N+(i-1)*M;
   s=x(1,i1:i2); %.*blackman(N)';
   
   wi=fft(s);
   %uwzgledniajac wlasnosci Matlabowej funkcji do fft:
   wi=wi/N*2;
   %stala normalizacji:
   wi=wi*lsb;
   wi=wi.*wi;
   %uwzglednienie kwantu f dla fft:
   wi=abs(wi)/df;
   %wartosc skuteczna sinusa to amplituda/sqrt(2) wiec:
   wi=wi/2;
   w=w+abs(wi);
end
w=w/steps;