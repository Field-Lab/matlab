function [Amplitudes]=NS_AmplitudesForPattern(DataPath,Channels,PatternNumber,MovieNumber,NS_GlobalConstants);
%The function gives back the amplitude value of stimulation current for
%each electrode.

FullName_status=[DataPath filesep 'status_m' num2str(MovieNumber)];
load(FullName_status);
FullName_pattern=[DataPath filesep 'pattern' num2str(PatternNumber) '_m' num2str(MovieNumber)];
load(FullName_pattern);
[patternTimes,traces] = plotStimulusTraces(Pattern, Status, Channels, [0 50], NS_GlobalConstants);
ST=size(traces)
for i=1:length(Channels)
    Trace=reshape(traces(1,i,:),1,ST(3));
    %figure(100);
    %plot(Trace);
    index=find(abs(Trace)==max(abs(Trace)));
    Amplitudes(i)=Trace(1,index(1));
end