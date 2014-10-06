N=40;
fs=44100;
f1=2000;
A=3;
fi=0;

t=[0:(N-1)]./fs;

s=A*sin(2*pi*f1*t+fi);
s0=round(s);

plot(t,s,t,s0)
