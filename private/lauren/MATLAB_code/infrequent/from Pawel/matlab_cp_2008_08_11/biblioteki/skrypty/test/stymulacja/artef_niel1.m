clear;

I=1e-6;
I0=35e-9;
Ce=6e-9;
Ut=52e-3;

a=I/Ce;
b=I0/Ce;
g=1/Ut;

dt=1e-6;
t=[0:100]*dt;

I=1e-6;
V0=0;
s=struct('I0',I0,'Ce',Ce,'Ut',Ut);
[y1,Ir1,Ic1]=artifact1niel(I,t,V0,s);
I=-I;
V0=y1(length(y1));
V0*1000
[y2,Ir2,Ic2]=artifact1niel(I,t,V0,s);

y=[y1 y2];
Ir=[Ir1 Ir2];
Ic=[Ic1 Ic2];
plot(y)

figure(1)
subplot(2,1,1)
plot(real(y))
grid on;
subplot(2,1,2)
plot(imag(y))
grid on

figure(2)
subplot(2,2,1)
plot(real(Ic1+Ir1))
grid on;
subplot(2,2,3)
plot(imag(Ic1+Ir1))
grid on
subplot(2,2,2)
plot(real(Ic2+Ir2))
grid on;
subplot(2,2,4)
plot(imag(Ic2+Ir2))
grid on

real(y(length(y)))*1000