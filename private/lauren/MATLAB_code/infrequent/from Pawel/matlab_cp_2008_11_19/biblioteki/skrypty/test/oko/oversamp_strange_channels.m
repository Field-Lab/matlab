%Oversampling - analiza "podejrzanych" kanalow, czyli:9,18,25,59,64.
cd I:\dane\2000\2000-12-12;
name='Data009';

nrchns=65;

samples=[s1 (s1+dl-1)];
fs=20000;

filtr_dl=20;
N=2;

marg_left=10; %tyle probek "na lewo" od poczatku spike'a zalicza sie do spike'a
marg_right=50;

%1. Kanal 9
channels=9;
s1=1000000;
s1=s1+1500000;
dl=1000000;
t=[1:dl]/fs;
figura=21;

read_param=struct('name',name,'header',206,'nrchns',nrchns,'channels',channels+1,'samples',samples);
s=readcnst(read_param);    
figure(figura);
subplot(2,2,1);
plot(t,s);
subplot(2,2,2);
f=[0:dl-1]/dl*fs;
plot(f,abs(fft(s)));
axis([0 10000 0 1000000]);

%1a). Oversampling:
N=2;
filtr_dl=16;
f_gr=0.98;
signal=s(1:2:dl); 

[y,filtr]=oversampling(signal,N,filtr_dl,f_gr);
y1=y(1,filtr_dl+1:filtr_dl+dl);
roznica=y1-s;

subplot(2,2,3);
plot(t,s);
axis([37.7 38.8 -600 600]);
subplot(2,2,4);
plot(t,roznica);
axis([37.7 38.8 -50 50]);

