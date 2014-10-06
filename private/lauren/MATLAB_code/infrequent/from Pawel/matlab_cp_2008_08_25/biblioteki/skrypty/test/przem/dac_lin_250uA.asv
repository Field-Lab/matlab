clear
cd H:/pliki/nauka/stymulacja/chip/testy/2006_03_15/
dane1=importdata('250uA_306.dat')';
dane2=importdata('250uA_1190.dat')';
dane3=importdata('250uA_4800.dat')';
dane4=importdata('250uA_19k.dat')';
%dane=importdata('1ua_76k_probabi.dat')';
figure(1);
clf;

R1=307;
R2=1190;
R3=4820;
R4=19170;
%R4=4684e3;

dac=[-127:127];
[slope1,lin_err1]=dac_lin3(dac,dane1);
[slope2,lin_err2]=dac_lin3(dac,dane2);
[slope3,lin_err3]=dac_lin3(dac,dane3);
[slope4,lin_err4]=dac_lin3(dac,dane4);

slope=slope1/R1*127
slope=slope2/R2*127
slope=slope3/R3*127
slope=slope4/R4*127

range=[-150 150 -1 1];

subplot(2,2,1);
plot(dac,lin_err1,'bd-');
axis(range);
a=text(0,1.1,['R=' num2str(R1) ',I=' num2str(slope1/R1*127*1e6) 'uA']);
set(a,'FontSize',16)
grid on;

subplot(2,2,2);
plot(dac,lin_err2,'bd-');
axis(range);
a=text(0,1.1,['R=' num2str(R2) ',I=' num2str(slope2/R2*127*1e6) 'uA']);
set(a,'FontSize',16)
grid on;

subplot(2,2,3);
plot(dac,lin_err3,'bd-');
axis(range);
a=text(0,1.1,['R=' num2str(R3) ',I=' num2str(slope3/R3*127*1e6) 'uA']);
set(a,'FontSize',16)
grid on;

subplot(2,2,4);
plot(dac,lin_err4,'bd-');
grid on;