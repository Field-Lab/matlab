clear;
figure(1)
R=[307 1190 4820 19170 74900 299e3 1190e3 4670e3 19200e3];
%*************************  1mA  **************************

cd H:/pliki/nauka/stymulacja/chip/testy/2006_03_15/
dane1=importdata('1mA_306.dat')';
dane2=importdata('1mA_1190.dat')';
dane3=importdata('1mA_4800.dat')';

R1=R(1,1);
R2=R(1,2);
R3=R(1,3);

dac=[-127:127];
[slope1,lin_err1]=dac_lin3(dac,dane1);
[slope2,lin_err2]=dac_lin3(dac,dane2);
[slope3,lin_err3]=dac_lin3(dac,dane3);

slope=slope1/R1*127
slope=slope2/R2*127
slope=slope3/R3*127

textx=-120;
texty=1.7;
textsize=14;
range=[-150 150 -1.5 1.5];

subplot(4,3,2);
plot(dac,lin_err1,'bd-');
axis(range);
a=text(textx,texty,['R=' num2str(R1) ', I=' num2str(slope1/R1*127*1e6) 'uA']);
set(a,'FontSize',textsize)
grid on;

subplot(4,3,3);
plot(dac,lin_err2,'bd-');
axis(range);
a=text(textx,texty,['R=' num2str(R2) ', I=' num2str(slope2/R2*127*1e6) 'uA']);
set(a,'FontSize',textsize)
grid on;

%*************************  250uA  **************************


%clear
cd H:/pliki/nauka/stymulacja/chip/testy/2006_03_15/
dane1=importdata('250uA_306.dat')';
dane2=importdata('250uA_1190.dat')';
dane3=importdata('250uA_4800.dat')';
%dane=importdata('1ua_76k_probabi.dat')';

R1=R(1,1);
R2=R(1,2);
R3=R(1,3);

dac=[-127:127];
[slope1,lin_err1]=dac_lin3(dac,dane1);
[slope2,lin_err2]=dac_lin3(dac,dane2);
[slope3,lin_err3]=dac_lin3(dac,dane3);

slope=slope1/R1*127
slope=slope2/R2*127
slope=slope3/R3*127

%textx=-100;
%texty=1.1;
textsize=14;
%range=[-150 150 -1 1];

subplot(4,3,4);
plot(dac,lin_err1,'bd-');
axis(range);
a=text(textx,texty,['R=' num2str(R1) ', I=' num2str(slope1/R1*127*1e6) 'uA']);
set(a,'FontSize',textsize)
grid on;

subplot(4,3,5);
plot(dac,lin_err2,'bd-');
axis(range);
a=text(textx,texty,['R=' num2str(R2) ', I=' num2str(slope2/R2*127*1e6) 'uA']);
set(a,'FontSize',textsize)
grid on;

subplot(4,3,6);
plot(dac,lin_err3,'bd-');
axis(range);
a=text(textx,texty,['R=' num2str(R3) ', I=' num2str(slope3/R3*127*1e6) 'uA']);
set(a,'FontSize',textsize)
grid on;


%*************************  60uA  **************************


%clear
cd H:/pliki/nauka/stymulacja/chip/testy/2006_03_13/
dane1=importdata('60ua_1200.dat')';
dane2=importdata('60ua_4800.dat')';
dane3=importdata('60ua_18k.dat')';
cd H:/pliki/nauka/stymulacja/chip/testy/2006_03_15/
dane4=importdata('60ua_75k.dat')';

R1=R(1,2);
R2=R(1,3);
R3=R(1,4);

dac=[-127:127];
[slope1,lin_err1]=dac_lin3(dac,dane1);
[slope2,lin_err2]=dac_lin3(dac,dane2);
[slope3,lin_err3]=dac_lin3(dac,dane3);

slope=slope1/R1*127
slope=slope2/R2*127
slope=slope3/R3*127

%textx=-100;
%texty=1.1;
textsize=14;
%range=[-150 150 -1 1];

subplot(4,3,7);
plot(dac,lin_err1,'bd-');
axis(range);
a=text(textx,texty,['R=' num2str(R1) ', I=' num2str(slope1/R1*127*1e6) 'uA']);
set(a,'FontSize',textsize)
grid on;

subplot(4,3,8);
plot(dac,lin_err2,'bd-');
plot(dac,lin_err2,'bd-');
axis(range);
a=text(textx,texty,['R=' num2str(R2) ', I=' num2str(slope2/R2*127*1e6) 'uA']);
set(a,'FontSize',textsize)
grid on;

subplot(4,3,9);
plot(dac,lin_err3,'bd-');
axis(range);
a=text(textx,texty,['R=' num2str(R3) ', I=' num2str(slope3/R3*127*1e6) 'uA']);
set(a,'FontSize',textsize)
grid on;


%*************************  16uA  **************************


%clear
cd H:/pliki/nauka/stymulacja/chip/testy/2006_03_10/
dane1=importdata('16ua_5k.dat')';
dane2=importdata('16ua_18k.dat')';
dane3=importdata('16ua_75k.dat')';
dane4=importdata('16ua_299k.dat')';
%dane=importdata('1ua_76k_probabi.dat')';

R1=4800;
R2=19100;
R3=74500;
R4=299000;
R1=R(1,3);
R2=R(1,4);
R3=R(1,5);
%R4=4684e3;

dac=[-127:127];
[slope1,lin_err1]=dac_lin3(dac,dane1);
[slope2,lin_err2]=dac_lin3(dac,dane2);
[slope3,lin_err3]=dac_lin3(dac,dane3);

slope=slope1/R1
slope=slope2/R2
slope=slope3/R3
%textx=-100;
%texty=0.7;
textsize=14;

%range=[-150 150 -1.4 0.6];

subplot(4,3,10);
plot(dac,lin_err1,'bd-');
axis(range);
a=text(textx,texty,['R=' num2str(R1) ', I=' num2str(slope1/R1*127*1e6) 'uA']);
set(a,'FontSize',textsize)
grid on;

subplot(4,3,11);
plot(dac,lin_err2,'bd-');
axis(range);
a=text(textx,texty,['R=' num2str(R2) ', I=' num2str(slope2/R2*127*1e6) 'uA']);
set(a,'FontSize',textsize)
grid on;

subplot(4,3,12);
plot(dac,lin_err3,'bd-');
axis(range);
a=text(textx,texty,['R=' num2str(R3) ', I=' num2str(slope3/R3*127*1e6) 'uA']);
set(a,'FontSize',textsize)
grid on;
