n=2000;
fs=20000;
t1=1;
t2=10*t1;
v1=-1;
v2=-0.1*v1;
marg=100;

n1=floor(fs*t1/1000);
n2=floor(fs*t2/1000);

waveform=zeros(1,n);
waveform(1,marg+1:marg+n1)=v1;
waveform(1,marg+n1+1:marg+n1+n2)=v2;

plot(waveform);

f=fft(waveform);
