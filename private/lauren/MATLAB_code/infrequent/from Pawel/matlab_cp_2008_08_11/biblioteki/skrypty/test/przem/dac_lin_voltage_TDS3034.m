clear
cd H:/pliki/nauka/stymulacja/chip/testy/2006_04_06/
cd H:/pliki/nauka/stymulacja/chip/testy/2006_04_21/

dane0=importdata('volt_r4_120k_b4_c1_dc500pos.dat')';
dane0=-importdata('volt_r3_120k_b5_c1_dc0neg.dat')';

dac=[-127:127];

dane1=[dane0(1,2:257)];
TDSoffset=[(dane0(1,1)+dane0(1,258))/2];

[slope1,lin_err1,dane]=dac_lin3(dac,dane1);

figure(1);
plot(dac,dane,'db-');
title('volt r6 120k b4 c1 dc500pos - data figure','FontSize',18,'FontWeight','demi');
xlabel('DAC (number)','FontSize',18,'FontWeight','demi');
ylabel('Voltage [V]','FontSize',18,'FontWeight','demi');
grid on;

dane2=dane1-TDSoffset;
[slope1,lin_err1,dane]=dac_lin3(dac,dane2);

h=gca;
set(h,'FontSize',16);
figure(2);
plot(dac,lin_err1,'bd-');
title('volt r6 120k b4 c1 dc500pos - linearity error','FontSize', 18,'FontWeight', 'demi');
xlabel('DAC (number)','FontSize', 18,'FontWeight', 'demi');
ylabel('linearity error','FontSize', 18,'FontWeight', 'demi');
h=gca;
set(h,'FontSize',16);

grid on;
