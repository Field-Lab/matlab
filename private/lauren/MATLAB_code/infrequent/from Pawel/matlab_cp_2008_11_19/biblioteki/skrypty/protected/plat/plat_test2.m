cd /home/pawel/pliki/dane/chip1_plat;

R=100000;

czytaj=1;

if czytaj==1
a0=importdata('plat0_chip1.dat');
a1=importdata('plat1_chip1.dat');
a2=importdata('plat2_chip1.dat');
a3=importdata('plat3_chip1.dat');
a4=importdata('plat4_chip1.dat');
a5=importdata('plat5_chip1.dat');
a8=importdata('plat8_chip1.dat');
a12=importdata('plat12_chip1.dat');
a57=importdata('plat57_chip1.dat');
a58=importdata('plat58_chip1.dat');
a59=importdata('plat59_chip1.dat');
a60=importdata('plat60_chip1.dat');
a61=importdata('plat61_chip1.dat');
a62=importdata('plat62_chip1.dat');
end

b(1,1:32)=mean(a0');
b(2,1:32)=mean(a1');
b(3,1:32)=mean(a2');
b(4,1:32)=mean(a3');
b(5,1:32)=mean(a4');
b(6,1:32)=mean(a5');
b(7,1:32)=mean(a8');
b(8,1:32)=mean(a12');
b(9,1:32)=mean(a57');
b(10,1:32)=mean(a58');
b(11,1:32)=mean(a59');
b(12,1:32)=mean(a60');
b(13,1:32)=mean(a61');
b(14,1:32)=mean(a62');

for i=1:14
    figure(1);
    subplot(3,5,i);
    plot(b(i,:)/R,'bd-');
    axis([0 31 -150e-8 0]);
    grid on;
    
    figure(2);
    subplot(3,5,i);
    plot(diff(b(i,:)/R),'bd-');
    axis([0 31 -8e-8 0]);
    grid on;
end
