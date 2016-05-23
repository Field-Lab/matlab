function bitstr=realtime_burst2(HoldWidth,HoldDelay);

FrameHeader=[1 0 1 0 1 1 1];
holdl = HoldDelay + HoldWidth;
Pause=zeros(1,holdl);
% finalnie: Frame=[FrameHeader Pause FrameData]

NumberOfFrames=20;
FrameLength=80;
%ChannelsOperated=[1];
NumberOfChannels=1;
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

channel=1; % channel 1 oznacza pierwszy kanal aktywny, niekoniecznie pierwszy kana? w og?le
delay=1; 
dac7b=20;
dac4b=1;
channel=1;
PulseShape=[Faza0' Faza1' Faza1' Faza2' Faza2' Faza3' Faza3' Faza4']';
realtime_sup_PH;

delay=11;
channel=1;
dac7b=20;
dac4b=1;
channel=1;
PulseShape=[Faza0' Faza1' Faza2' Faza3' Faza4']';
realtime_sup_PH;

bitstr=RejectSpaces(num2str(DataStream));
