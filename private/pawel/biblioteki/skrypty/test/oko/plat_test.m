cd /home/pawel/pliki/dane/;

a61=importdata('plat61_chip1.dat');
a62=importdata('plat62_chip1.dat');
a61_100=importdata('plat61_chip1_100.dat');

b61=mean(a61');
b62=mean(a62');
b61_100=mean(a61_100');

t=[0:31];

figure(1);
plot(t,b61,'bd-',t,b62,'rd-',t,b61_100,'gd-');
%plot(t,b61_100,'bdfigure(2);
figure(2);
%plot(t,b61_100);

c61=std(a61(4,:))
c61_100=std(a61_100(4,:))

d61=diff(b61);
d61_100=diff(b61_100);
std(d61)
std(d61_100)