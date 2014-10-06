function [StimChannels,Amplitudes]=NS_StimulatedChannels(DataPath,PatternNumber,MovieNumber,Channels,NS_GlobalConstants);
%WARNING!! This is function that works only for typical scan, when ALL the
%active channels have the same range. This function was built specifically
%to work with scan data recorded with incorrect information of the current
%ranges - the erro is with the chip address: for Labview it is from 11111
%down to 11000 (chips 1 to 8) and for Matlab it was always from 11000 to
%11111. THIS HAS TO BE CHANGED (PH, 2010-12-06)
%The output value is the maximum absolut value of amplitude over all the
%Channels

FullName_status=[DataPath filesep 'status_m' num2str(MovieNumber)];
load(FullName_status);

r=0;
for i=1:length(Status.ChannelsStatus)
    a=Status.ChannelsStatus(i).range;
    if a>r
        r=a;
    end
end

for i=1:length(Status.ChannelsStatus)
    Status.ChannelsStatus(i).range=r;
    %Status.ChannelsStatus(i).range
end    

FullName_pattern=[DataPath filesep 'pattern' num2str(PatternNumber) '_m' num2str(MovieNumber)];
load(FullName_pattern);
[patternTimes,traces] = plotStimulusTraces(Pattern,Status,Channels,[0 40],NS_GlobalConstants);
Straces=size(traces);
Traces=zeros(Straces(2),Straces(3));
Traces(:,:)=traces(1,:,:);
a=max(abs(Traces')); % max value for each channel
StimChannels=find(a>0);
Amplitudes=a(StimChannels);