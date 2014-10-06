Ampl=1;
T=1000;
tphase=400;
dt=1;
N=floor(T/dt);
t=[dt:dt:T];

delay1=100;
tau=100;

i1=zeros(1,N);
i2=i1;

i1(201:300)=-0.8;
i1(301:400)=0.8;

i2(101:200)=0.6;
i2(201:300)=-1;
i2(301:400)=0.4;

fs=18;
x1=0;
x2=1;
y1=-1.2;
y2=1.2;
%subplot(1,2,1);
figura=1;
figure(figura);
a=plot(t/1000,i1);
%set(a,'Color',[51/256,153/256,102/256]);
set(a,'LineWidth',2)
axis([x1 x2 y1 y2]);
grid on;
h=gca;
set(h,'FontSize',fs)
set(h,'XTick',[0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1]);
xlabel('t [ms]');
ylabel('I_{stym} [\muA]');
print(figura,'-dpng','IstimBiphasic');

Ampl=1;
for j=1:7                    
        %t_disch=200e-6;
        A1=(0.4+0.05*j)*Ampl
        A2=-Ampl;
        A3=-A1-A2;        
        i2(101:200)=A1;
        i2(201:300)=A2;
        i2(301:400)=A3;
        w1(j,:)=i2;
        size(w1);
end

figura=2;
figure(figura);
y1=-1.2;
y2=1.2;
a=plot(t/1000,w1(1,:),t/1000,w1(2,:),t/1000,w1(3,:),t/1000,w1(4,:),t/1000,w1(5,:),t/1000,w1(6,:),t/1000,w1(7,:));
%a=plot(t/1000,i2)
%set(a,'Color',[51/256,153/256,102/256]);
set(a,'LineWidth',2)
axis([x1 x2 y1 y2]);
grid on
h=gca;
set(h,'FontSize',fs)
set(h,'XTick',[0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1]);
xlabel('t [ms]');
ylabel('I_{stym} [\muA]');

a=legend('A1=0.45','A1=0.50','A1=0.55','A1=0.60','A1=0.65','A1=0.70','A1=0.75');
print(figura,'-dpng','IstimTriphasic');
