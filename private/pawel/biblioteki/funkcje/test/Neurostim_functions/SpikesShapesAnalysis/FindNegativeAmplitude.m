function NegAmplitude=FindNegativeAmplitude(Waveform);
%Traces: MxN (M - number of channels, N - number of smaples)
[x1,x2]=min(Waveform,[],2);
NegAmplitude=x1;

%threshold=x1/2;
%find(Waveform-threshold<0)
%HalfAmplitudeWidth=length(find(Waveform-threshold<0));