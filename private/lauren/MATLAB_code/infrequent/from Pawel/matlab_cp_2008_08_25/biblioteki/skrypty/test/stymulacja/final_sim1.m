cd H:\pliki\nauka\stymulacja\symulacje\2005_06_20;

a=importdata('5u_sym_nodisch_el.out');
s0=(a(:,2)-2.5)*1000;
a=importdata('5u_sym_100us_el.out');
s1=(a(:,2)-2.5)*1000;
a=importdata('5u_sym_200us_el.out');
s2=(a(:,2)-2.5)*1000;
a=importdata('5u_sym_300us_el.out');
s3=(a(:,2)-2.5)*1000;

a=importdata('5u_nonsym_nodisch_el.out');
z0=(a(:,2)-2.5)*1000;
a=importdata('5u_nonsym_100us_el.out');
z1=(a(:,2)-2.5)*1000;
a=importdata('5u_nonsym_200us_el.out');
z2=(a(:,2)-2.5)*1000;
a=importdata('5u_nonsym_300us_el.out');
z3=(a(:,2)-2.5)*1000;

l0=ones(1,1501)*1.2;
l1=ones(1,1501)*(-1.2);

t=[-100:1400];
h=plot(t,s0,t,s1,t,s2,t,s3,t,z0,t,z1,t,z2,t,z3,t,l0,t,l1);
col2=[0 0.5 0];
set(h(1),'LineWidth',2)
set(h(1),'Color','Blue')
set(h(2),'LineWidth',2)
set(h(2),'Color','Blue')
set(h(3),'LineWidth',2)
set(h(3),'Color','Blue')
set(h(4),'LineWidth',2)
set(h(4),'Color','Blue')
set(h(5),'LineWidth',2)
set(h(5),'Color',col2)
set(h(6),'LineWidth',2)
set(h(6),'Color',col2)
set(h(7),'LineWidth',2)
set(h(7),'Color',col2)
set(h(8),'LineWidth',2)
set(h(8),'Color',col2)
set(h(9),'LineWidth',2)
set(h(9),'Color','Red')
set(h(9),'LineStyle','--')
set(h(10),'LineWidth',2)
set(h(10),'Color','Red')
set(h(10),'LineStyle','--')

grid on
axis([-100 1400 -5 5]);
xlabel('time [us]');
ylabel('chip input [mV]')

h=gca;
set(h,'FontSize',16)