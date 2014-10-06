function Amplitude=NS_StimulationAmplitude(DataPath,PatternNumber,MovieNumber,Channels,NS_GlobalConstants);
%The output value is the maximum absolut value of amplitude over all the
%Channels

FullName_status=[DataPath filesep 'status_m' num2str(MovieNumber)];
load(FullName_status);
FullName_pattern=[DataPath filesep 'pattern' num2str(PatternNumber) '_m' num2str(MovieNumber)];
load(FullName_pattern);
[patternTimes,traces] = plotStimulusTraces(Pattern,Status,Channels,[0 40],NS_GlobalConstants);
Amplitude=max(max(traces));