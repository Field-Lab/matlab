ActiveBitsPerFrame=13*NumberOfChannelsOperated;
ClkMask=zeros(1,FrameLength);
ClkMask(HoldLength+FrameHeaderLength+1:HoldLength+FrameHeaderLength+ActiveBitsPerFrame)=2; %2 means real+time data

Frame=zeros(1,FrameLength);
Frame(HoldLength+FrameHeaderLength+1:13:HoldLength+FrameHeaderLength+1+13*(NumberOfChannelsOperated-1))=1;
StreamLength=20; % in sampling periods
DataStream=zeros(1,StreamLength*FrameLength);
ClkStream=zeros(1,StreamLength*FrameLength);
TriggerConnectStream=zeros(1,StreamLength*FrameLength);
TriggerRecordStream=zeros(1,StreamLength*FrameLength);

for i=1:StreamLength
    DataStream((i-1)*FrameLength+1:i*FrameLength)=Frame;
    ClkStream((i-1)*FrameLength+1:i*FrameLength)=ClkMask;
    TriggerConnectStream((i-1)*FrameLength+TriggerConnectTiming:(i-1)*FrameLength+TriggerConnectTiming+TrigerDuration-1)=1;
    TriggerRecordStream((i-1)*FrameLength+TriggerRecordTiming:(i-1)*FrameLength+TriggerRecordTiming+TrigerDuration-1)=1;
end

Faza0=[0 0 0 1 2]; % record connect polarity dac6b dac4b
Faza1=[0 1 0 1 2];
Faza2=[0 1 1 1 3];
Faza3=[0 1 0 1 1];
Faza4=[0 0 0 1 1];

amplitude2=1;
channel=1;
amplitude=20;
delay=2;
PulseShape=[Faza0' Faza1' Faza1' Faza2' Faza2' Faza3' Faza3' Faza4']';
Cadence_sim_realtime_sup;

amplitude2=3;
channel=1;
amplitude=40;
delay=12;
PulseShape=[Faza0' Faza1' Faza2' Faza3' Faza4']';
Cadence_sim_realtime_sup;

amplitude2=1;
channel=3;
amplitude=40;
delay=5;
PulseShape=[Faza0' Faza1' Faza1' Faza2' Faza2' Faza3' Faza3' Faza4']';
Cadence_sim_realtime_sup;

amplitude2=2;
channel=3;
amplitude=60;
delay=15;
PulseShape=[Faza0' Faza1' Faza2' Faza3' Faza4']';
Cadence_sim_realtime_sup;

channel=8;
amplitude=60;
delay=10;
amplitude2=3;
PulseShape=[Faza0' Faza1' Faza2' Faza3' Faza4']';
%Cadence_sim_realtime_sup;

FullDataStream=[FullDataStream DataStream];
FullClkStream=[FullClkStream ClkStream];
FullTriggerConnectStream=[FullTriggerConnectStream TriggerConnectStream];
FullTriggerRecordStream=[FullTriggerRecordStream TriggerRecordStream];