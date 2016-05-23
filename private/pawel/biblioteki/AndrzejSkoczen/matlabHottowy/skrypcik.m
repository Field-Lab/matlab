ActiveChannels=[1:4]; % liczone od 1
HW=40;
HD=0;
FrameLength=120;

s3=realtime_burst1(HW,HD);
Channel=1; % drugi kanal ze wszystkich aktywnych!
offset=HW+HD+7+(Channel-1)*14;
ramka=5; %liczone od 1

word=s3([offset+1:offset+14]+(ramka-1)*FrameLength)