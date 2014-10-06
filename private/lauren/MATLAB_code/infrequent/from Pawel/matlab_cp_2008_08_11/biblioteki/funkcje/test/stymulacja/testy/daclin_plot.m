function w=daclin_plot(range,

dane1=importdata('1ua_19k.dat')';
dane2=importdata('1ua_76k.dat')';
dane3=importdata('1ua_299k.dat')';
dane4=importdata('1ua_1200k.dat')';
%dane=importdata('1ua_76k_probabi.dat')';

R1=19100;
R2=74500;
R3=299000;
R4=1189e3;
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

figure(3)

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