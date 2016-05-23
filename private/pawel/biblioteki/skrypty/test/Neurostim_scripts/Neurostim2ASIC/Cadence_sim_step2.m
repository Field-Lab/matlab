for p=1:length(PulseShape)
    Phase=PulseShape(p,:)
    PhaseData=[Phase(1:3) de2bi(Phase(4)*amplitude,6,'left-msb') de2bi(Phase(5)*amplitude2,4,'left-msb')]
    
    FrameNumber=p+delay;
    DataStream((FrameNumber-1)*FrameLength+HoldLength+(channel-1)*13+1:(FrameNumber-1)*FrameLength+HoldLength+channel*13)=PhaseData;
end