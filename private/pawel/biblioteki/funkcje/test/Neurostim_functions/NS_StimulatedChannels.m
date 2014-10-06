function [StimChannels,Amplitudes]=NS_StimulatedChannels(DataPath,PatternNumber,MovieNumber,Channels,NS_GlobalConstants);
%The output value is the maximum absolut value of amplitude over all the
%Channels

FullName_status=[DataPath filesep 'status_m' num2str(MovieNumber)];
load(FullName_status);
Status
%for i=1:512
%    Status.ChannelsStatus(i).range=Status.ChannelsStatus(12).range
%end
FullName_pattern=[DataPath filesep 'pattern' num2str(PatternNumber) '_m' num2str(MovieNumber)];
load(FullName_pattern);
[patternTimes,traces] = plotStimulusTraces(Pattern,Status,Channels,[0 40],NS_GlobalConstants);
Straces=size(traces);
Traces=zeros(Straces(2),Straces(3));
Traces(:,:)=traces(1,:,:);
a=max(abs(Traces')); % max value for each channel
StimChannels=find(a>0);
Amplitudes=a(StimChannels);