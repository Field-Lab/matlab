cd H:\pliki\nauka\stymulacja\chip\testy\2006_05_24;
clear;
fs=14;

a1=importdata('volt_r5_299k_b5_c1_ground.dat')';
a2=importdata('volt_r5_75k_b5_c1_ground.dat')';
a3=importdata('volt_r5_51k_b5_c1_ground.dat')';

cd H:\pliki\nauka\stymulacja\chip\testy\2006_05_23;
a4=importdata('volt_r5_19100_b5_c1_ground.dat')';

dac=[-127:127];

dane1=[a1(1,1:127) (a1(1,128)+a1(1,129))/2 a1(1,130:256)];
dane2=[a2(1,1:127) (a2(1,128)+a2(1,129))/2 a2(1,130:256)];
dane3=[a3(1,1:127) (a3(1,128)+a3(1,129))/2 a3(1,130:256)];
dane4=[a4(1,1:127) (a4(1,128)+a4(1,129))/2 a4(1,130:256)];

a=plot(dac,dane1,'bd-',dac,dane3,'kd-',dac,dane4,'rd-')
grid on
axis([-140 140 -2 2])

h=gca;
set(h,'FontSize',fs);
set(h,'LineWidth',2)
a=xlabel('DAC setting');
set(a,'FontSize',fs);
a=ylabel('output voltage [V]');
set(a,'FontSize',fs);
a-legend('RL=300k','RL=50k','RL=20k');