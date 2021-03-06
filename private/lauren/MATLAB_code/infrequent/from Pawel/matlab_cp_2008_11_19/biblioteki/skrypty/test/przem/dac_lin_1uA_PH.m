clear
cd H:/pliki/nauka/stymulacja/chip/testy/2006_03_16/
dane1=importdata('1ua_75k_b4_c1.dat')';
dane2=importdata('1ua_299k_b4_c1.dat')';
dane3=importdata('1ua_1200k_b4_c1.dat')';
dane4=importdata('1ua_1200k_b4_c1.dat')';

figure(1);
clf;

R1=74500;
R2=299000;
R3=1189e3;
R4=4684e3;
%R4=4684e3;

dac=[-127:127];

fit=1;
if fit==0
    [slope1,lin_err1,dane1,offset1]=dac_lin5(dac,dane1);
    [slope2,lin_err2,dane2,offset2]=dac_lin5(dac,dane2);
    [slope3,lin_err3,dane3,offset3]=dac_lin5(dac,dane3);
    [slope4,lin_err4,dane4,offset4]=dac_lin5(dac,dane4);
else
    [slope1,lin_err1,dane1,offset1]=dac_lin5_abs_sum(dac,dane1);
    [slope2,lin_err2,dane2,offset2]=dac_lin5_abs_sum(dac,dane2);
    [slope3,lin_err3,dane3,offset3]=dac_lin5_abs_sum(dac,dane3);
    [slope4,lin_err4,dane4,offset4]=dac_lin5_abs_sum(dac,dane4);
end

slope=slope1/R1
slope=slope2/R2
slope=slope3/R3
slope=slope4/R4

osie1=[-140 140 -1.2 1.2];
osie2=[-140 140 -1.2 1.2];
dane=dane3;
lin_err=lin_err3;
R=R3;
fs=18;

%subplot(3,2,1)
figure(1)
plot(dac,dane/R*1e6,'kd-');
axis(osie1)
grid on;
h=gca;
set(h,'FontSize',fs);
set(h,'LineWidth',2)
a=xlabel('DAC value');
set(a,'FontSize',fs);
a=ylabel('output current [\muA]');
set(a,'FontSize',fs);
%subplot(3,2,2)
figure(5)
plot(dac,lin_err,'kd-');
axis(osie2)
grid on;
h=gca;
set(h,'FontSize',fs);
set(h,'LineWidth',2)
a=xlabel('DAC value');
set(a,'FontSize',fs);
a=ylabel('linearity error [LSB]');
set(a,'FontSize',fs);

figure(3)
plot(dane4)
