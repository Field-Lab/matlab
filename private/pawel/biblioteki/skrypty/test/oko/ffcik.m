clc;
clear;
f=6000;
A=1;
n=2000;
f_s=[15000 14000 13000 12000 11000 10000];
for i=1:length(f_s)
fs=f_s(i)
t=0:n-1;
t=t*1/fs;
x=A*cos(2*pi*f*t);
wabs=abs(fft(x));
wa=wabs(1:(n/2)+1)/n*2;
p=(n/2)+1
df=fs/n;
freq=0:(p-1);
freq=freq*df;
subplot(2,3,i);
plot(freq,wa);
end 