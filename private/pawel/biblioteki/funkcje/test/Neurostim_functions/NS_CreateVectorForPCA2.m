function OutputTraces=NS_CreateVectorForPCA2(Traces,TracesChannels,Pattern,NS_GlobalConstants);
%Pattern is simply a subarray of the original "patterns_out" array, given
%by "ReadPAtternDataChunk function. (PH, 2008-09-15)
%Traces - 3-dimensional array: traces*electrodes*samples.
%TracesChannels - includes the number of electrodes to which data in
%"Traces" array correspond. The length of this array must be identical to
%second dimension of the "Traces" array.
%Pulse - array of the size 5*N, where N is the length of the pulse in
%sampling periods, and the 5 rows are:
%row 1 - the DAC values;
%row 2 - values of "record" signal;
%row 3 - values of "connect" signal;
%row 4 - values of "discharge" signal;
%row 5 - values of "hold" signal.
%Status - describes the status of the chips and individual channels.
%FileName;
%Channel;
%PDChunkIndex;

ChipAddresses=NS_GlobalConstants.ChipAddresses;
NumberOfChannelsPerChip=NS_GlobalConstants.NumberOfChannelsPerChip;
CurrentRanges=NS_GlobalConstants.CurrentRanges;
Fs=NS_GlobalConstants.SamplingFrequency;

STraces=size(Traces);
MaskInit=zeros(1,STraces(3));
OutputTraces=[];
%Pulse=zeros(length(Channels,

%Find stimulated channels
Pattern
StimChannels=[];
for i=1:length(Pattern)
    StimChannels=[StimChannels Pattern(i).channel];
end

%StimChannels=Pattern(:).channel;

for i=1:length(TracesChannels)
    Channel=TracesChannels(i);
    clear FindChannel;
    FindChannel=find(StimChannels==Channel)
    if length(FindChannel)==1
        data=Traces(:,i,:);
        Pulse=Pattern(FindChannel).data;
        Disconnect=Pulse(2,:);
        Mask=MaskInit;
        Mask(1:length(Disconnect))=Disconnect;
        WhichSamples=find(Mask==0);
        signal=reshape(Traces(:,i,WhichSamples),STraces(1),length(WhichSamples));
    elseif length(FindChannel)==0
        signal=reshape(Traces(:,i,:),STraces(1),STraces(3));
    else
        error('channel found more than once in this pattern')
    end
    OutputTraces=[OutputTraces signal];
end