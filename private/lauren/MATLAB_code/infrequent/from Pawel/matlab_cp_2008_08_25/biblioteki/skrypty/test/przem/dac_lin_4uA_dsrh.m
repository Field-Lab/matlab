clear
cd H:/pliki/nauka/stymulacja/chip/testy/2006_03_10/
dane1=importdata('4ua_18k_b4_c1.dat')';
dane2=importdata('4ua_75k_b4_c1.dat')';
dane3=importdata('4ua_299k_b4_c1.dat')';
dane4=importdata('4ua_299k_b4_c1.dat')';
%dane=importdata('1ua_76k_probabi.dat')';
figure(1);
clf;

R1=19100;
R2=74500;
R3=299000;
R4=1189e3;
%R4=4684e3;

dac=[-127:127];
[slope1,lin_err1,dane11]=dac_lin3(dac,dane1);
[slope2,lin_err2]=dac_lin3(dac,dane2);
[slope3,lin_err3]=dac_lin3(dac,dane3);
[slope4,lin_err4]=dac_lin3(dac,dane4);

slope=slope1/R1
slope=slope2/R2
slope=slope3/R3
slope=slope4/R4

subplot(2,2,1);
plot(dac,dane11,'bd-');
grid on;

subplot(2,2,2);
plot(dac,lin_err2,'bd-');
grid on;

subplot(2,2,3);
plot(dac,lin_err3,'bd-');
grid on;
subplot(2,2,4);
plot(dac,lin_err4,'bd-');
grid on;

figure(2)
[a,b,c]=plotyy(dac,dane11/R1,dac,dane11)
x=get(a(2),'Ylim')
set(a(1),'Ylim',x/R1)
grid on
set(b,'LineWidth',3)
set(c,'LineWidth',3)
%set(b,'LineStyle','')
set(get(a(1),'Ylabel'),'String','Left Y-axis')
set(get(a(2),'Ylabel'),'String','Right Y-axis')