clear
cd H:/pliki/nauka/stymulacja/chip/testy/2006_03_13/

dane1=importdata('250nA_299k_b4_c1.dat')';
dane2=importdata('250nA_1200k_b4_c1.dat')';
dane3=importdata('250nA_4700k_b4_c1.dat')';
dane4=importdata('250nA_19200k_b4_c1.dat')';

figure(1);
clf;

R1=299000;
R2=1189e3;
R3=4684e3;
R4=19200e3;
%R4=4684e3;

dac=[-127:127];
[slope1,lin_err1]=dac_lin3(dac,dane1);
[slope2,lin_err2]=dac_lin3(dac,dane2);
[slope3,lin_err3]=dac_lin3(dac,dane3);
[slope4,lin_err4]=dac_lin3(dac,dane4);

slope=slope1/R1
slope=slope2/R2
slope=slope3/R3
slope=slope4/R4

subplot(2,2,1);
plot(dac,lin_err1,'bd-');
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