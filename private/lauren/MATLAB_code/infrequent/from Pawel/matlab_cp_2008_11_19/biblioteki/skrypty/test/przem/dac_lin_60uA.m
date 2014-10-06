clear
cd H:/pliki/nauka/stymulacja/chip/testy/2006_03_13/
dane1=importdata('60ua_1200.dat')';
dane2=importdata('60ua_4800.dat')';
dane3=importdata('60ua_18k.dat')';
cd H:/pliki/nauka/stymulacja/chip/testy/2006_03_15/
dane4=importdata('60ua_75k.dat')';
%dane4=importdata('025ua_19200k.dat')';
%dane=importdata('1ua_76k_probabi.dat')';
figure(1);
clf;

R1=1200;
R2=4800;
R3=19100;
R4=74500;

dac=[-127:127];
[slope1,lin_err1]=dac_lin3(dac,dane1);
[slope2,lin_err2]=dac_lin3(dac,dane2);
[slope3,lin_err3]=dac_lin3(dac,dane3);
[slope4,lin_err4]=dac_lin3(dac,dane4);

slope=slope1/R1*127
slope=slope2/R2*127
slope=slope3/R3*127
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