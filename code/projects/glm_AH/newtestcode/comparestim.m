figure; 
frames = 400:640;
msec   = frames*8-400*8;
plot(msec,squeeze(testmoviecone.matrix(20,20,frames)),'r'); hold on; plot(msec,squeeze(testmoviekernel.matrix(20,20,frames+1)),'b')
xlabel('msecs');
ylabel('unitless stim');
title('Red: Full Cone Model, Blue: Convolved Time Kernel');

orient landscape
eval(sprintf('print -dpdf FullConeStim_vs_TimeKernelStim.pdf'))