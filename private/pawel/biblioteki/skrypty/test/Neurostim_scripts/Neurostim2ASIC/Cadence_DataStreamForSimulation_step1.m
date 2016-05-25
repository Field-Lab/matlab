NumberOfChannels=64; %numbered from 0 to N-1!

SamplingPeriod=25e-6;

FrameLength=1000;
ChannelsOperated=ones(1,NumberOfChannels);
NumberOfOperatedChannels=length(find(ChannelsOperated==1))
%NeutralFrame:
Frame=zeros(1,FrameLength);
Frame(HoldLength+1:13:HoldLength+1+13*(NumberOfOperatedChannels-1))=1;
StreamLength=10; % in sampling periods
DataStream=zeros(1,StreamLength*FrameLength);
for i=1:StreamLength
    DataStream((i-1)*FrameLength+1:i*FrameLength)=Frame;
end
plot(DataStream,'bd-')
grid on