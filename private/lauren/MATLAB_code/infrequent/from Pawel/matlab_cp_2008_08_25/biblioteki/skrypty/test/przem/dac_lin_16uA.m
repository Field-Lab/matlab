clear
cd H:/pliki/nauka/stymulacja/chip/testy/2006_03_10/
dane1=importdata('16ua_5k.dat')';
dane2=importdata('16ua_18k.dat')';
dane3=importdata('16ua_75k.dat')';
dane4=importdata('16ua_299k.dat')';
%dane=importdata('1ua_76k_probabi.dat')';
figure(1);
clf;

R1=4800;
R2=19100;
R3=74500;
R4=299000;
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
textx=-100;
texty=0.7;
textsize=14;

range=[-150 150 -1.4 0.6];
subplot(4,3,1);
plot(dac,lin_err1,'bd-');
axis(range);
a=text(textx,texty,['R=' num2str(R1) ', I=' num2str(slope1/R1*127*1e6) 'uA']);
set(a,'FontSize',textsize)
grid on;

subplot(4,3,2);
plot(dac,lin_err2,'bd-');
axis(range);
a=text(textx,texty,['R=' num2str(R2) ', I=' num2str(slope2/R2*127*1e6) 'uA']);
set(a,'FontSize',textsize)
grid on;

subplot(4,3,3);
plot(dac,lin_err3,'bd-');
axis(range);
a=text(textx,texty,['R=' num2str(R3) ', I=' num2str(slope3/R3*127*1e6) 'uA']);
set(a,'FontSize',textsize)
grid on;
%subplot(2,2,4);
%plot(dac,lin_err4,'bd-');
%grid on;