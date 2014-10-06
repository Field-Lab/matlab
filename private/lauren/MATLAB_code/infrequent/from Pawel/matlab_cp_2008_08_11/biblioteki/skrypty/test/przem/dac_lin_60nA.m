clear
cd H:/pliki/nauka/stymulacja/chip/testy/2006_03_13/
dane1=importdata('60na_1200k.dat')';
dane2=importdata('60na_4700k.dat')';
dane3=importdata('60na_19200k.dat')';
%dane4=importdata('025ua_19200k.dat')';
%dane=importdata('1ua_76k_probabi.dat')';
figure(1);
clf;

R1=1189e3;
R2=4684e3;
R3=19200;
%R4=4684e3;

dac=[-127:127];
[slope1,lin_err1]=dac_lin3(dac,dane1);
[slope2,lin_err2]=dac_lin3(dac,dane2);
[slope3,lin_err3]=dac_lin3(dac,dane3);
%[slope4,lin_err4]=dac_lin3(dac,dane4);

slope=slope1/R1*127
slope=slope2/R2*127
slope=slope3/R3*127
%slope=slope4/R4

subplot(2,2,1);
plot(dac,lin_err1,'bd-');
grid on;

subplot(2,2,2);
plot(dac,lin_err2,'bd-');
grid on;

subplot(2,2,3);
plot(dac,lin_err3,'bd-');
grid on;

%subplot(2,2,4);
%plot(dac,lin_err4,'bd-');
%grid on;