Fs = 20000;
L = 400000;    
NFFT = 2^nextpow2(L); % Next power of 2 from length of y
Y = fft(s2,NFFT)/L;
f = Fs/2*linspace(0,1,NFFT/2+1);

% Plot single-sided amplitude spectrum.

figure
plot(f,2*abs(Y(1:NFFT/2+1)))
axis ([0 100 0 0.00002])
title('Single-Sided Amplitude Spectrum of y(t)')
xlabel('Frequency (Hz)')
ylabel('|Y(f)|')

