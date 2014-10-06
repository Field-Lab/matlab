cd /home2/pawel/oko/2000-12-12;

nrchns=65;
s1=239700;
dl=4000;
samples=[s1 (s1+dl-1)];
channels=17;
%channels=31;
name='Data009';
read_param=struct('name',name,'header',206,'nrchns',nrchns,'channels',channels,'samples',samples);
s=readcnst(read_param);
time=[1:dl]/20;
figure(1);
plot(s);
%axis([0 time(dl) -1100 700]);
grid on;

figure(2);
channels=65;
s1=150000;
samples=[s1 (s1+dl-1)];
read_param=struct('name',name,'header',206,'nrchns',nrchns,'channels',channels,'samples',samples);
s65=readcnst(read_param);
subplot(3,1,1);
plot(time,s65);
axis([0 time(dl) -1100 700]);
grid on;

channels=31;
s1=298600;
samples=[s1 (s1+dl-1)];
read_param=struct('name',name,'header',206,'nrchns',nrchns,'channels',channels,'samples',samples);
s31=readcnst(read_param);
subplot(3,1,2);
plot(time,s31);
axis([0 time(dl) -1100 700]);
grid on;

channels=17;
s1=239700;
samples=[s1 (s1+dl-1)];
read_param=struct('name',name,'header',206,'nrchns',nrchns,'channels',channels,'samples',samples);
s17=readcnst(read_param);
subplot(3,1,3);
plot(time,s17);
axis([0 time(dl) -1100 700]);
grid on;
