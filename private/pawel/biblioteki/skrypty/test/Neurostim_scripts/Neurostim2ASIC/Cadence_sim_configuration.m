ChannelsOperated=ones(1,NumberOfChannels);
NumberOfChannelsOperated=length(find(ChannelsOperated==1));

StimAnalog=ones(1,NumberOfChannelsOperated);
StimDigital=ones(1,NumberOfChannelsOperated);
RecAnalog=ones(1,NumberOfChannelsOperated);
RecMode=ones(1,NumberOfChannelsOperated);
GroundInput=ones(1,NumberOfChannelsOperated);

Frame=zeros(1,FrameLength);

Frame(HoldLength+FrameHeaderLength+1:5:HoldLength+FrameHeaderLength+1+5*(NumberOfChannelsOperated-1))=StimAnalog;
Frame(HoldLength+FrameHeaderLength+2:5:HoldLength+FrameHeaderLength+2+5*(NumberOfChannelsOperated-1))=StimDigital;
Frame(HoldLength+FrameHeaderLength+3:5:HoldLength+FrameHeaderLength+3+5*(NumberOfChannelsOperated-1))=RecAnalog;
Frame(HoldLength+FrameHeaderLength+4:5:HoldLength+FrameHeaderLength+4+5*(NumberOfChannelsOperated-1))=RecMode;
Frame(HoldLength+FrameHeaderLength+5:5:HoldLength+FrameHeaderLength+5+5*(NumberOfChannelsOperated-1))=GroundInput;

Stream=[Stream Frame];