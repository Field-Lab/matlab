function EI=NS_CalculateEI(Traces);
%Calculates the average traces for each channel.
%Input:
%Traces - array of size NxCxT, where:
%N - number of traces for each channel;
%C - number of channels;
%T - length of single trace.
STraces=size(Traces);
EI=zeros(STraces(2),STraces(3));

for i=1:STraces(2) %for each channel...
    ChannelData=Traces(:,i,:);   
    ChannelDataBis=reshape(ChannelData,STraces(1),STraces(3));
    SData=size(ChannelDataBis);
    if SData(1)~=1
        Waveform=mean(ChannelDataBis);
    else
        Waveform=ChannelDataBis;
    end
    EI(i,:)=Waveform;    
end