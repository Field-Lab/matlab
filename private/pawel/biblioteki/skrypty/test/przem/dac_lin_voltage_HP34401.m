clear
cd H:/pliki/nauka/stymulacja/chip/testy/2006_04_06/

dane0=importdata('volt_r4_120k_b4_c1_dc100pos.dat')';
dac=[-127:127];

dane1=[dane0(1,2:257)];
TDSoffset=[(dane0(1,1)+dane0(1,258))/2];

[slope1,lin_err1,dane]=dac_lin3(dac,dane1);

figure;
plot(dac,dane,'-db');
title('volt r4 120k b4 c1 dc0pos - data figure');
xlabel('DAC (number)');
ylabel('Voltage [V]');
grid on;

dane2=dane1-TDSoffset;
[slope1,lin_err1,dane]=dac_lin3(dac,dane2);

figure;
plot(dac,lin_err1,'bd-');
title('volt r4 120k b4 c1 dc0pos - linearity error','FontSize', 12,'FontWeight', 'demi');
xlabel('DAC (number)','FontSize', 12,'FontWeight', 'demi');
ylabel('linearity error','FontSize', 12,'FontWeight', 'demi');

grid on;


