ChannelIndex=channel
for p=1:length(PulseShape)
    Phase=PulseShape(p,:);
    PhaseData=[Phase(1:3) de2bi(Phase(4)*dac7b,7,'left-msb') de2bi(Phase(5)*dac4b,4,'left-msb')];
    
    FrameNumber=p+delay
    Index=HoldDelay+HoldWidth+length(FrameHeader)+(FrameNumber-1)*FrameLength+(ChannelIndex-1)*14
    DataStream(Index+1:Index+14)=PhaseData;    
end