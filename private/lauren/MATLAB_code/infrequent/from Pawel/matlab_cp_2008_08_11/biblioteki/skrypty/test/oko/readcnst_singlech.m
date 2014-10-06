cd /home2/pawel/oko/2000-12-12;

nrchns=65;
s1=3239700+134190;
dl=320;
samples=[s1 (s1+dl-1)];
channels=65;
%channels=31;
name='Data009';
fs=20000;
read_param=struct('name',name,'header',206,'nrchns',nrchns,'channels',channels,'samples',samples);
s=readcnst(read_param);
time=[1:dl]/20;


%nrchns=65;
%channels=65;
%dl=500;
%s1=149950;
%samples=[s1 (s1+dl-1)];
%read_param=struct('name',name,'header',206,'nrchns',nrchns,'channels',channels,'samples',samples);
%s1=readcnst(read_param);
%time=[1:dl]/20000*1000;


figure(2);
plot(time,s);
axis([0 16 -1200 700]);

xlabel('time [miliseconds]');
ylabel('signal value');
text(12,560,'2000/12/12');
text(12,495,'Data009');
text(12,430,'channel 64');
grid on;

figure(3);
f=fft(s)/dl*2*dl/90;
cz=[0:dl-1]/dl*fs;
plot(cz,abs(f));
axis([0 fs/2 0 300]);
text(8100,240,'2000/12/12');
text(8100,225,'Data009');
text(8100,210,'channel 64');
xlabel('frequency [Hz]');
ylabel('fourier coefficient');
grid on;
