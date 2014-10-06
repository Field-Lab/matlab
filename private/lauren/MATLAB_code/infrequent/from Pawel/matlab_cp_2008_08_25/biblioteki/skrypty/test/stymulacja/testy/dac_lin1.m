cd H:/pliki/nauka/stymulacja/chip/testy/2006_03_01/
a=importdata('PROBA.dat')';
b=importdata('PROBA2.dat')';

dac=[-127:127];
dane=[b(1,128:-1:2) (a(1,1)+b(1,1))/2 a(1,2:1:128)];

[p,s] = POLYFIT(dac,dane,1);
y = POLYVAL(p,dac);
figure(1)
plot(dac,dane);
grid on
figure(2)
plot(dac,y-dane);
grid on