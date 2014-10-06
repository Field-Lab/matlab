function [x] = minuscz (t,x,f0)
fs = 20000;             %czêstotliwoœæ próbkowania
% f0 - czêstotliwoœæ notch
fn = fs/2;              %czêstotliwoœæ Nyquista
freqRatio = f0/fn;      
notchWidth = 0.01;       %szerokoœæ notch

%Liczenie zer
zeros = [exp( sqrt(-1)*pi*freqRatio ), exp( -sqrt(-1)*pi*freqRatio )];

%Liczenie biegunów
poles = (1-notchWidth) * zeros;

figure;
zplane(zeros.', poles.');

b = poly( zeros ); %# Get moving average filter coefficients
a = poly( poles ); %# Get autoregressive filter coefficients

figure;
freqz(b,a,2000,fs)

%filtracja sygna³u x
y = filter(b,a,x);
figure
% pierwszy wykres
subplot(2,1,1)
plot(t,x);
xlabel('Sec')
axis tight
title('Oryginalny sygna³')
%drugi wykres
subplot(2,1,2)
plot(t,y);
xlabel('kHz')
title('Sygna³ przefiltrowany')