rawDataDir = '/Volumes/Acquisition/Data/2014-12-10-SB1test/';
datarun = 'data001'; 
rawDataPath = [rawDataDir datarun];
rawFile = edu.ucsc.neurobiology.vision.io.RawDataFile(rawDataPath);
header = rawFile.getHeader();
totalSamples = header.getNumberOfSamples(); 
channels = 1:(header.getNumberOfElectrodes-1); %-1 because first channel has TTL pulses
NS_GlobalConstants = NS_GenerateGlobalConstants(length(channels));

%%
plottingOn = 1;
Fs = 20000;                   % Sampling frequency
T = 1/Fs;                     % Sample time
L = 10000;                    % Length of signal
t = (0:L-1)*T;                % Time vector
NFFT = 2^nextpow2(L);         % Next power of 2 from length of y
f = Fs/2*linspace(0,1,NFFT/2+1);

idx = 0; 
% Reading in the data 

startSample = 10000;
nSamplesToRead = totalSamples-startSample;
fftData = zeros(totalSamples/10000-1,163,512); 
if plottingOn; figure; end
while nSamplesToRead>0
    idx = idx + 1; 

    rawData = double(rawFile.getData(startSample,min(10000,nSamplesToRead)));
    nSamplesToRead = nSamplesToRead - 10000;
    startSample = startSample + 10000;

    y = rawData(:,2:end) - repmat(mean(rawData(:,2:end)),10000,1);
    Y = fft(y,NFFT)/L;

   
    fftData(idx,:,:) = 2*abs(Y(3:165,:)); 
    fprintf('\n%d samples to read; %0.1f%% complete\n',nSamplesToRead,100*startSample/totalSamples)
    if plottingOn
%         subplot(2,1,1); cla;
%         plot(1000*t(1:100),rawData(1:100,2:end))
%         title('raw signal from all channels')
%         xlabel('time (milliseconds)');
%         title(sprintf('%d samples to read',nSamplesToRead));
%         % Plot single-sided amplitude spectrum.
%         subplot(2,1,2); cla;
        plot(f(3:165),2*abs(Y(3:165,:)))
        title('Single-Sided Amplitude Spectrum of y(t)')
        xlabel('Frequency (Hz)');
        ylabel('|Y(f)|');
        xlim([3 200]); 
        pause(0.05)
    end
end
rawFile.close; 
%% 
Y = squeeze(mean(fftData));
figure; 
plot(f(3:165),Y);
title(sprintf('Single-Sided Amplitude Spectrum for %s',datarun))
xlabel('Frequency (Hz)');
ylabel('|Y(f)|');
xlim([3 200]);
grid on; 
freqRange = f(3:165); 
channelAverage = mean(Y,2);
noiseFrequency = freqRange(find(channelAverage == max(channelAverage)));

% figure; 
% for chan = 1:512
%     cla;
%     % Plot single-sided amplitude spectrum.
%     plot(f(3:165),Y(:,chan));
%     title(sprintf('Single-Sided Amplitude Spectrum of raw data on channel %d',chan))
%     xlabel('Frequency (Hz)');
%     ylabel('|Y(f)|'); 
%     xlim([3 200]); 
%     pause(0.1);
% end

