Ampl=1;
T=1000;
tphase=400;
dt=1;
N=floor(T/dt);
t=[dt:dt:T];

delay1=100;
tau=100;

v1=zeros(1,N);
v2=v1;
v3=v1;
start=delay1+1;
stop=N;

for i=start:stop
    v1(i)=-Ampl*(exp(-(i-start+1)/(tau/dt))-exp(-(stop-start+1)/(tau/dt)));    
end

start=delay1+tphase+1;
stop=N;
for i=start:stop
    v2(i)=Ampl*(exp(-(i-start+1)/(tau/dt))-exp(-(stop-start+1)/(tau/dt)));
end

start=delay1+1;
stop=delay1+tphase;

for i=start:stop
    v3(i)=-Ampl;   
end

x1=-0.1;
x2=1.1;
y1=-1.2;
y2=0.2;
subplot(4,1,2);
a=plot(t/1000,v3);
set(a,'Color',[51/256,153/256,102/256]);
set(a,'LineWidth',2)
axis([x1 x2 y1 y2]);
a=gca;
set(a,'Visible','off');

subplot(2,1,2)
x1=-0.1;
x2=1.1;
y1=-1.2;
y2=1.2;
a=plot(t/1000,v1+v2)
set(a,'Color',[51/256,153/256,102/256]);
set(a,'LineWidth',2)
axis([x1 x2 y1 y2]);
a=gca;
set(a,'Visible','off');
%a=gca;
%h=gca;
%fs=20;
%set(h,'FontSize',fs);
%set(h,'XTick',[1e1 1e2 1e3 1e4 1e5 1e6]);
%set(h,'LineWidth',2);
%xlabel('czas [ms]');
%ylabel('I_{stim} [k\Omega]');