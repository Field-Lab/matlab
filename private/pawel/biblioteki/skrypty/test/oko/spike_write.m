readcnst_channel;

fs=20000;
spike=s65(1,1:1800);
t=[1:length(spike)]/fs;
w=trapez_blackman(length(spike),50);
spike=spike.*w/2000000;
y=[t' spike']';

cd /home/pawel;
fid = fopen('spike.txt','w');
fprintf(fid,'%12.8f  %12.8f\n',y);
fclose(fid);
