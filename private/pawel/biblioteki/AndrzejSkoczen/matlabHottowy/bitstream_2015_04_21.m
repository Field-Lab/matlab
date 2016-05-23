%fid = fopen('C:\pawel\nauka\Neurostim-3\datastream2.txt','w');
HW = 40;
HD = 0;
nrch = 64;

Masks=zeros(1,nrch);
ActiveChannels=[1:4];
DL1=length(ActiveChannels)*14;
Masks1=Masks; Masks1(ActiveChannels)=1;

ActiveChannels=[61:64];
DL2=length(ActiveChannels)*14;
Masks2=Masks; Masks2(ActiveChannels)=1;

ActiveChannels=[1];
DL3=length(ActiveChannels)*14;
Masks3=Masks; Masks3(ActiveChannels)=1;

ActiveChannels=[64];
DL4=length(ActiveChannels)*14;
Masks4=Masks; Masks4(ActiveChannels)=1;

bitstr = [SoftReset, '__',HoldTiming(HW,HD),'__',TriggerDelay(10,5,6),'__'];
bitstr = [bitstr,ChannelMaskBlockPH(DL1,Masks1),'__',ChannelConfig(nrch,nrch*14),'__'];
bitstr = [bitstr,ChannelMaskBlockPH(DL1,Masks1),'__',realtime_burst1(HW,HD),'__'];
bitstr = [bitstr,ChannelMaskBlockPH(DL2,Masks2),'__',realtime_burst1(HW,HD),'__'];
bitstr = [bitstr,ChannelMaskBlockPH(DL3,Masks),'__',realtime_burst2(HW,HD),'__'];
bitstr = [bitstr,ChannelMaskBlockPH(DL3,Masks),'__',realtime_burst2(HW,HD),'__'];
%fprintf(fid,'%s',bitstr);
%fclose(fid);