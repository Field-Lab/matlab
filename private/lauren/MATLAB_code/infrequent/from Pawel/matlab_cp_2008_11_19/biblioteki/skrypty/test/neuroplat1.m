cd /home/pawel/pliki/nauka

a=importdata('pasmo1.dat');
b=mean(a);
c=importdata('czest.dat');
c1=c(1,2:length(c));
b1=b(1,2:28);

f=[20 30 40 50 70 120 320 720 920 1200 1500 1800 2000 2400 3000 4000 5000];
g=[350 580 620 700 710 710 710 710 700 700 680 640 600 545 480 380 320];

b1=b1*max(max(g))/max(max(b1));

loglog(c1,b1,'bd-',f,g,'r*-');
grid on;
axis([10 10000 100 1000]);
legend('internal calibration','sinus generator');
xlabel('frequency');
ylabel('gain (normalized)');
