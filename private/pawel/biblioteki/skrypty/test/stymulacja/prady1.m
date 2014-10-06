t0=0.1e-3;
t1=1e-4; %czas trwania pierwszej czesci impulsu
t2=1e-4; %t2 - czas trwanai drugiej czesci
t3=1e-4;
dt=1e-6;

A1=0;
A2=-0.8;
A3=+0.8;
B1=0.6;
B2=-1;
B3=0.4;

prad=zeros(1,10000);
A=prad;
B=prad;

for i=1001:2000
    A(i)=A1;
    B(i)=B1;
end
    
for i=2001:3000
    A(i)=A2;
    B(i)=B2;
end

for i=3001:4000
    A(i)=A3;
    B(i)=B3;
end
t=[1:10000]/10000;

figure(1)
a=plot(t,A);
set(a,'LineWidth',2);
axis([0 1 -1.2 1.2]);
grid on;
h=gca;
xlabel('czas [ms]');
ylabel('I_{stym} [\muA]');
set(h,'FontSize',20)
set(h,'LineWidth',2);
set(h,'XTick',[0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1]);

figure(2)
a=plot(t,B);
set(a,'LineWidth',2);
axis([0 1 -1.2 1.2]);
grid on;
h=gca;
xlabel('czas [ms]');
ylabel('I_{stym} [\muA]');
set(h,'FontSize',20)
set(h,'LineWidth',2);
set(h,'XTick',[0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1]);
%t=[dt:dt:T]*1000;