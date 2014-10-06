% cd /mnt/data4/dane/2000/2000-12-12; %linux
cd I:\dane\2000\2000-12-12; %windows
clear;

s1=150000;
dl=4000;
spike_location=29.6; %in miliseconds
spike_location=56.45; %in miliseconds
spike_location=9.3;
spike_start=6;
spike_stop=14;
%spike_location=175.6;

sp_loc=spike_location;

nrchns=65;
%s1=3239700+134190;
%dl=320;
samples=[s1 (s1+dl-1)];
channels=65;
%channels=31;
name='Data009';
fs=20000;
read_param=struct('name',name,'header',206,'nrchns',nrchns,'channels',channels,'samples',samples);
s=readcnst(read_param);
time=[1:dl]/20;

cd I:\dane\2004\12maja\ForPawel\2003-08-06-0-016;
s=importdata('electrode167.txt')';

%Nadprobkowanie:
fp=20;
t0=[0:length(s)-1]/fp;
czest_gr=0.95;

filtr_dl=36;
N=5;
[y,filtr]=oversampling(s,N,filtr_dl,czest_gr);
y0=y(1,1+filtr_dl/2*N:length(y)-filtr_dl/2*N+N-1);
clear y;
t1=[0:length(y0)-1]/(N*fp);

'sdfsdgsg'

%break;
signal=s(1:2:dl);
t2=[0:length(signal)-1]/fp*2;
N=10;
[y1,filtr]=oversampling(signal,N,filtr_dl,czest_gr);
y2=y1(1,1+filtr_dl/2*N:length(y1)-filtr_dl/2*N+N-1);
clear y1;
t3=[0:length(y2)-1]/(N*fp)*2;
'fwf'
%break;
figure(2)
plot(t0,s,'b-',t1,y0,'g-',t2,signal,'r-',t3,y2,'k-');
legend('raw 20kHz','20kHz->100kHz','undersampled 10kHz','10kHz->100kHz');
grid on;