clear
cd H:/pliki/nauka/stymulacja/chip/testy/2006_03_15/
dane1=importdata('1mA_306_b4_c1.dat')';
dane2=importdata('1mA_1190_b4_c1.dat')';
dane3=importdata('1mA_4800_b4_c1.dat')';

%dane=importdata('1ua_76k_probabi.dat')';
figure(1);
clf;

R1=306;
R2=1190;
R3=4800;

%R4=4684e3;

dac=[-127:127];
[slope1,lin_err1]=dac_lin3(dac,dane1);
[slope2,lin_err2]=dac_lin3(dac,dane2);
[slope3,lin_err3]=dac_lin3(dac,dane3);

slope=slope1/R1*12700
slope=slope2/R2*12700
slope=slope3/R3*12700

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
grid on;