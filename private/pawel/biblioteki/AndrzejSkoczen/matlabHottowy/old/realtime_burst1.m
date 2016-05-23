function bitstr = realtime_burst1(HoldDelay,HoldWidth)
FrameHeader=[1 0 1 0 1 1 1];
Pause=zeros(1,HoldDelay+HoldWidth);
% finalnie: Frame=[FrameHeader Pause FrameData]
NumberOfFrames=20;
FrameLength=120;
NumberOfChannels=4;
ActiveBitsPerFrame=14*NumberOfChannels; %Operated;
NeutralFrameData=zeros(1,ActiveBitsPerFrame);

DataStream=zeros(1,NumberOfFrames*FrameLength);
for i=1:NumberOfFrames
    DataStream((i-1)*FrameLength+1:(i-1)*FrameLength+7) = FrameHeader;
end

Faza0=[0 0 0 1 2]; % record connect polarity dac7b dac4b
Faza1=[0 1 0 1 2];
Faza2=[0 1 1 1 3];
Faza3=[0 1 0 1 1];
Faza4=[0 0 0 1 1];

%kolejne impulsy:
bitstr = '';
channel=1;
delay=1;
dac7b=20;
dac4b=1;
PulseShape=[Faza0' Faza1' Faza1' Faza2' Faza2' Faza3' Faza3' Faza4']';
bitstr =[bitstr,RejectSpaces(realtime_sup(DataStream,channel,delay,dac7b,dac4b,PulseShape,HoldDelay,HoldWidth,FrameHeader,FrameLength))];

channel=2;
delay=3;
dac7b=20;
dac4b=1;
PulseShape=[Faza0' Faza1' Faza1' Faza2' Faza2' Faza3' Faza3' Faza4']';
bitstr =[bitstr,RejectSpaces(realtime_sup(DataStream,channel,delay,dac7b,dac4b,PulseShape,HoldDelay,HoldWidth,FrameHeader,FrameLength))];

channel=3;
delay=5;
dac7b=20;
dac4b=1;
PulseShape=[Faza0' Faza1' Faza1' Faza2' Faza2' Faza3' Faza3' Faza4']';
bitstr =[bitstr,RejectSpaces(realtime_sup(DataStream,channel,delay,dac7b,dac4b,PulseShape,HoldDelay,HoldWidth,FrameHeader,FrameLength))];

channel=4;
delay=7;
dac7b=20;
dac4b=1;
PulseShape=[Faza0' Faza1' Faza1' Faza2' Faza2' Faza3' Faza3' Faza4']';
bitstr =[bitstr,RejectSpaces(realtime_sup(DataStream,channel,delay,dac7b,dac4b,PulseShape,HoldDelay,HoldWidth,FrameHeader,FrameLength))];

channel=1;
delay=11;
dac7b=10;
dac4b=1;
PulseShape=[Faza0' Faza1' Faza1' Faza2' Faza2' Faza3' Faza3' Faza4']';
bitstr =[bitstr,RejectSpaces(realtime_sup(DataStream,channel,delay,dac7b,dac4b,PulseShape,HoldDelay,HoldWidth,FrameHeader,FrameLength))];

channel=2;
delay=13;
dac7b=20;
dac4b=2;
PulseShape=[Faza0' Faza1' Faza1' Faza2' Faza2' Faza3' Faza3' Faza4']';
bitstr =[bitstr,RejectSpaces(realtime_sup(DataStream,channel,delay,dac7b,dac4b,PulseShape,HoldDelay,HoldWidth,FrameHeader,FrameLength))];
%disp(bitstr);
end
