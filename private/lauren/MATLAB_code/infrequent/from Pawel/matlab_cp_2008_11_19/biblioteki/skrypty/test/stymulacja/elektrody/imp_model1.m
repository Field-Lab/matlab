c=12e+7;
f=[20:1:20000];

Rs=35000;
Re=2e7;

beta=0.76;
zcpe=c*(i*2*pi*f).^(-beta);
a=zcpe+Re;

a1=zcpe*Re;
z1=a1./a;
z=z1+Rs;
pomiar(1,:)=[20 30 40 60 90 130 270 520 1020 2020 4020 7000 10000 20000 40000 80000 150000];
pomiar(2,:)=[2850 2120 1700 1320 970 735 427 263 165 106 72 59 52 45 42 38 38];
pomiar(3,:)=-[59.4 60.1 60 64 64.5 63.8 64.2 62.5 58 52 43 34 30 22 19 19 25];

figura=10;
figure(figura);
subplot(2,1,1);
loglog(f,abs(z)/1000,pomiar(1,:),pomiar(2,:),'g*');
grid on;
subplot(2,1,2);
semilogx(f,angle(z)*57,pomiar(1,:),pomiar(3,:),'g*')
grid on;

%beta=1;
%model=struct('Ce',3e7,'beta',0.64,'Re',2.57e6,'Rs',5130);
model=struct('Ce',c,'beta',beta,'Re',Re,'Rs',20000);
t0=0.1e-3;
t1=1e-4; %czas trwania pierwszej czesci impulsu
A1=0.6e-6; %amplituda pierwszej czesci
%A1=0;
t2=1e-4; %t2 - czas trwanai drugiej czesci
A2=-1e-6; %A2 - amplituda drugiej czesci
t3=1e-4;
A3=-A1-A2;
%A3=A3*0.95;
t_delay=1e-5;
t_disch=2e-4;
%t_disch=0;
%A2=0;
T=1e-3; %T - czas symulacji
dt=0.1e-6; %dt - krok czasowy
t=[dt:dt:T]*1000;

w1=artifact2(model,0,t0,t1,A1,t2,A2,t3,A3,t_delay,t_disch,T,dt);
w2=artifact2(model,0,t0,t1,A1,t2,A2,t3,A3,t_delay,0,T,dt);
%w2=artifact1(model,t0,t1,A1,t2,A2,T,1e-7);
figure(7)
plot(t,w1*2,t,w2*2);
%axis([0 2 -0.002 0.002]);
%plot([1:2000],w1,[1:2000],w2(10:10:20000));
grid on;
figure(3);
plot(t,w1-w2);
grid on;