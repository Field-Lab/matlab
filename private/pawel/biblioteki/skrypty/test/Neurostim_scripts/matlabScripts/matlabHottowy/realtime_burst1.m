function bitstr=realtime_burst1(HoldWidth,HoldDelay);

FrameHeader=[1 0 1 0 1 1 1];
holdl = HoldDelay + HoldWidth;
Pause=zeros(1,holdl);
% finalnie: Frame=[FrameHeader Pause FrameData]

NumberOfFrames=20;
FrameLength=120;
%ChannelsOperated=[1:4];
NumberOfChannels=4;
ActiveBitsPerFrame=14*NumberOfChannels;
NeutralFrameData=zeros(1,ActiveBitsPerFrame);

DataStream=zeros(1,NumberOfFrames*FrameLength);
for i=1:NumberOfFrames
    DataStream((i-1)*FrameLength+1:(i-1)*FrameLength+7+holdl) = [FrameHeader Pause];
end

Faza0=[0 0 0 1 2]; % record connect polarity dac7b dac4b
Faza1=[0 1 0 1 2];
Faza2=[0 1 1 1 3];
Faza3=[0 1 0 1 1];
Faza4=[0 0 0 1 1];

%kolejne impulsy:

channel=1;
delay=1;
dac7b=20;
dac4b=1;
channel=1;
PulseShape=[Faza0' Faza1' Faza1' Faza2' Faza2' Faza3' Faza3' Faza4']';
realtime_sup_PH;

channel=2;
delay=3;
dac7b=20;
dac4b=1;
channel=1;
PulseShape=[Faza0' Faza1' Faza1' Faza2' Faza2' Faza3' Faza3' Faza4']';
realtime_sup_PH;

channel=3;
delay=5;
dac7b=20;
dac4b=1;
channel=1;
PulseShape=[Faza0' Faza1' Faza1' Faza2' Faza2' Faza3' Faza3' Faza4']';
realtime_sup_PH;

channel=4;
delay=7;
dac7b=20;
dac4b=1;
channel=1;
PulseShape=[Faza0' Faza1' Faza1' Faza2' Faza2' Faza3' Faza3' Faza4']';
realtime_sup_PH;

channel=1;
delay=11;
dac7b=10;
dac4b=1;
channel=1;
PulseShape=[Faza0' Faza1' Faza1' Faza2' Faza2' Faza3' Faza3' Faza4']';
realtime_sup_PH;

channel=2;
delay=12;
dac7b=20;
dac4b=2;
channel=1;
PulseShape=[Faza0' Faza1' Faza1' Faza2' Faza2' Faza3' Faza3' Faza4']';
realtime_sup_PH;

bitstr = RejectSpaces(num2str(DataStream));
