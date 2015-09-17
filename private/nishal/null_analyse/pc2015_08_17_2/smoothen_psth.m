function PSTH=smoothen_psth(PSTH)
fftPSTH = fft(PSTH);

freq = 2*pi*((1:length(fftPSTH))-1)/length(fftPSTH);
filterfft=zeros(length(fftPSTH),1)';
x=0.99;
filterfft(freq<pi*(1-x)  | freq >pi*(1+x))=1;
fftPSTHnew=fftPSTH.*filterfft;
PSTHnew = ifft(fftPSTHnew);

% figure;
% plot(PSTH,'r');
% hold on;
% plot(PSTHnew,'b');


PSTH=PSTHnew;

end