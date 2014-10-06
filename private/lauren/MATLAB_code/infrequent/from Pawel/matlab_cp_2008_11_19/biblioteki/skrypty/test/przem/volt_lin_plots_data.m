cd H:/pliki/nauka/stymulacja/chip/testy/2006_04_24/

a=importdata('volt_r0_299k_b5_c1_ground.dat');
b=importdata('volt_r1_299k_b5_c1_ground.dat');
c=importdata('volt_r2_299k_b5_c1_ground.dat');
d=importdata('volt_r3_299k_b5_c1_ground.dat');
e=importdata('volt_r4_299k_b5_c1_ground.dat');
f=importdata('volt_r5_299k_b5_c1_ground.dat');
g=importdata('volt_r6_299k_b5_c1_ground.dat');
h=importdata('volt_r7_299k_b5_c1_ground.dat');

figure(1);

subplot(2,2,1);
plot(a,'bd-');
grid on;

subplot(2,2,2);
plot(b,'bd-');
grid on;

subplot(2,2,3);
plot(c,'bd-');
grid on;

subplot(2,2,4);
plot(d,'bd-');
grid on;