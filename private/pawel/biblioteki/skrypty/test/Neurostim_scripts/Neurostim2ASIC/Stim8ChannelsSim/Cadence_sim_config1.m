FrameLength=80;
ChannelsOperated=[1:NumberOfChannels];
NumberOfChannelsOperated=length(ChannelsOperated);

ActiveBitsPerFrame=NumberOfChannelsOperated*5;
ClkMask=zeros(1,FrameLength);
ClkMask(HoldLength+FrameHeaderLength+1:HoldLength+FrameHeaderLength+NumberOfChannels*5)=1; % 1 means configuration

StimAnalog=zeros(1,NumberOfChannels);
StimAnalog(1,ChannelsOperated)=1;
StimDigital=zeros(1,NumberOfChannels);
StimDigital(1,ChannelsOperated)=1;
RecAnalog=ones(1,NumberOfChannels);
RecMode=ones(1,NumberOfChannels);
GroundInput=ones(1,NumberOfChannels);

Frame=zeros(1,FrameLength);

Frame(HoldLength+FrameHeaderLength+1:5:HoldLength+FrameHeaderLength+1+5*(NumberOfChannels-1))=StimAnalog;
Frame(HoldLength+FrameHeaderLength+2:5:HoldLength+FrameHeaderLength+2+5*(NumberOfChannels-1))=StimDigital;
Frame(HoldLength+FrameHeaderLength+3:5:HoldLength+FrameHeaderLength+3+5*(NumberOfChannels-1))=RecAnalog;
Frame(HoldLength+FrameHeaderLength+4:5:HoldLength+FrameHeaderLength+4+5*(NumberOfChannels-1))=RecMode;
Frame(HoldLength+FrameHeaderLength+5:5:HoldLength+FrameHeaderLength+5+5*(NumberOfChannels-1))=GroundInput;

FullDataStream=[FullDataStream Frame];
FullClkStream=[FullClkStream ClkMask];

TriggerConnectTiming=65;
TriggerRecordTiming=75;

FullTriggerConnectStream=[FullTriggerConnectStream zeros(1,length(Frame))];
FullTriggerRecordStream=[FullTriggerRecordStream zeros(1,length(Frame))];