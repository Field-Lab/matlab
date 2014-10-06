n=2200;
fp=2000;
t=[1:n]/fp;

f1=340;
f2=340*345/335

a1=1;
a2=0.3;

s=a1*sin(2*pi*f1*t)+a2*sin(2*pi*f2*t);
plot(t,s)
grid on;