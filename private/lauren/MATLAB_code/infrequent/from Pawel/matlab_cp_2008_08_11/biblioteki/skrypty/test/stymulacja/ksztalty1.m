Ampl=1;
T=800e-6;
tphase=100e-6;
dt=1e-6;
delay1=100e-6;
delay2=100e-6;
tau=25e-6;

I=[dt:dt:T];

I1=I;
I1(floor(delay1/dt+1):floor(delay1/dt+tphase/dt))=-Ampl;
I1(floor(delay1/dt+tphase/dt+1):floor(delay1/dt+2*tphase/dt))=Ampl;
%plot(I1);

I2=I;
I2(floor(delay1/dt+1):floor(delay1/dt+tphase/dt))=-Ampl;
I2(floor(delay1/dt+tphase/dt+delay2/dt+1):floor(delay1/dt+delay2/dt+2*tphase/dt))=Ampl;
plot(I2);

I3=I;
I3(floor(delay1/dt+1):floor(delay1/dt+tphase/dt))=-Ampl;
I3(floor(delay1/dt+tphase/dt+1):floor(delay1/dt+6*tphase/dt))=Ampl/5;
plot(I3);

I4=I;
I4(floor(delay1/dt+1):floor(delay1/dt+tphase/dt))=-Ampl;
start=floor(delay1/dt+tphase/dt+1)
stop=floor(delay1/dt+2*tphase/dt);
for i=start:stop
  I4(i)=4.487*Ampl*(exp(-(i-start+1)/(tau/dt))-exp(-(stop-start+1)/(tau/dt)));
end
mean(I4)
%I4(floor(delay1/dt+tphase/dt+1):floor(delay1/dt+2*tphase/dt))=exp(-;
plot(I4);

x1=-100;
x2=900;
y1=-1.2*Ampl;
y2=-y1;

subplot(3,2,1);
a=plot(I1);
set(a,'Color',[51/256,153/256,102/256]);
set(a,'LineWidth',2)
axis([x1 x2 y1 y2]);
a=gca;
set(a,'Visible','off');

subplot(3,2,3);
a=plot(I2);
set(a,'Color',[51/256,153/256,102/256]);
set(a,'LineWidth',2)
axis([x1 x2 y1 y2]);
a=gca;
set(a,'Visible','off');

subplot(3,2,5);
a=plot(I3);
set(a,'Color',[51/256,153/256,102/256]);
set(a,'LineWidth',2)
axis([x1 x2 y1 y2]);
a=gca;
set(a,'Visible','off');

subplot(1,2,2);
a=plot(I4);
set(a,'Color',[51/256,153/256,102/256]);
set(a,'LineWidth',2)
axis([x1 x2 3.5*y1 3.5*y2]);
a=gca;
set(a,'Visible','off');