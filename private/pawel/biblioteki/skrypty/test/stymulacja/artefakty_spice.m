%cd /home/pawel/pliki/stymulacja/luty2005;
cd H:\pliki\nauka\stymulacja\symulacje\2005_02_24\

a1=importdata('artifact_bal_disch_100u.out');
a2=importdata('artifact_bal_disch_200u.out');
a3=importdata('artifact_bal_disch_300u.out');
a4=importdata('artifact_bal_disch_400u.out');
a5=importdata('artifact_bal_disch_500u.out');

b1=importdata('artifact_bal_no_disch.out');
c1=importdata('artifact_unbal.out');

t=b1(:,1)*1000;

art1=a1(:,2)*1000;
art2=a2(:,2)*1000;
art3=a3(:,2)*1000;
art4=a4(:,2)*1000;
art5=a5(:,2)*1000;

brt1=b1(:,2)*1000;
crt1=c1(:,2)*1000;
figure(1);
ap=plot(t,brt1,'k-.',t,crt1,'k--',t,art1,'k-',t,art2,'k-',t,art3,'k-',t,art4,'k-')
set(ap(1,1),'LineWidth',2);
set(ap(2,1),'LineWidth',2);
set(ap(3,1),'LineWidth',2);
set(ap(4,1),'LineWidth',2);
set(ap(5,1),'LineWidth',2);
set(ap(6,1),'LineWidth',2);

axis([0 2 -10 10]);
grid on;
h=gca;
set(h,'FontSize',16);
%legend('symmetric, no discharging','non-symmetric','symmetric, discharged');
a=xlabel('time [ms]');
ylabel('input artifact [mV]');
grid on;

curr=importdata('current_short.out');
a=importdata('artifact_bal_no_disch_short.out');
art0=a(:,2)*1000;
t=a(:,1)*1000;

figure(2)

subplot(2,1,1);
ap=plot(t,curr(:,2)*1000000,'k-')
set(ap,'LineWidth',2);
h=gca;
set(h,'FontSize',16);
a=xlabel('time [ms]');
ylabel('stimulus current [uA]');
axis([0 0.4 -6 6]);
grid on;

subplot(2,1,2);
ap=plot(t,art0,'k-')
set(ap,'LineWidth',2);
h=gca;
set(h,'FontSize',16);
a=xlabel('time [ms]');
ylabel('time course for Va [mV]');
axis([0 0.4 -500 1000]);
grid on;
