fid = fopen('C:\home\Pawel\nauka\granty\Symfonia2013\realizacja\chip\datastream.txt','w');
HW = 40;
HD = 4;
nrch = 64;
CM0=zeros(1,64);

bitstr = [SoftReset, '__',HoldTiming(HW,HD),'__',TriggerDelay(10,5,6),'__'];

DL=56;
CM=CM0;
CM(61:64)=1;
bitstr = [bitstr,ChannelMask_PH(nrch,DL,CM),'__',ChannelConfig(nrch,DL),'__'];
bitstr = [bitstr,ChannelMask_PH(nrch,DL,CM),'__',realtime_burst1(HW,HD),'__'];

DL=14;
CM=CM0;
CM(64)=1;
bitstr = [bitstr,ChannelMask_PH(nrch,DL,CM),'__',realtime_burst2(HW,HD),'__'];

DL=14;
CM=CM0;
CM(63)=1;
bitstr = [bitstr,ChannelMask_PH(nrch,DL,CM),'__',realtime_burst2(HW,HD),'__'];

%AnalogSettings(DACs.Gain,50)]
fprintf(fid,'%s',bitstr);
fclose(fid);

%HoldTiming(HW=40,HD=4)
% TriggerDelay(TD1=10, TD2=5, TD2=5)
% ChannelConfig, CHCL=56, CHC=<1111.....1111>
% ChannelMask, DL=56, CM=<0000.....00001111>
% RealTimeFrame*20: m??j skrypt Cadence_sim_realtime__v2_burst1

% ChannelMask, DL=14, CM=<0000.....0001>
% RealTimeFrame*20: m??j skrypt Cadence_sim_realtime__v2_burst2

% ChannelMask, DL=14, CM=<0000.....0010>
% RealTimeFrame*20: m??j skrypt Cadence_sim_realtime__v2_burst2

