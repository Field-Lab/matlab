function HalfAmplitudeWidth=FindNegativeHalfAmplitudeWidth(Waveform);
%Traces: MxN (M - number of channels, N - number of smaples)
N=10;
WaveformUpsampled=Upsample(Waveform,N);
[x1,x2]=min(WaveformUpsampled,[],2);
threshold=x1/2;
%find(Waveform-threshold<0)
HalfAmplitudeWidth=length(find(WaveformUpsampled-threshold<0))/N;