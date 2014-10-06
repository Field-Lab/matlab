cd /home/pawel/pliki/dane/
a4M=importdata('plat0_chip1_0_minus5V_4M.dat');
cd /home/pawel/pliki/dane/chip2_plat;
a220k=importdata('plat0_chip2_220kOhm.dat');
a3M=importdata('plat0_chip2_3_7M.dat');
values=[0:31];

a4M=a4M/3700000;
a220k=a220k/220000;
a3M=a3M/3700000;

a1=mean(a220k')*1e6;
a2=mean(a3M')*1e6;
a3=mean(a4M')*1e6;

numer1=1;
numer2=2;

figure(numer1);
figure(numer2);

clf(numer1);
clf(numer2);

[p,s]=polyfit(values,a1,1);
f=polyval(p,values);
LSB=abs(p(1))
s=abs(std(a1-f)/LSB)

[p,s]=polyfit(values(1:20),a2(1:20),1);
f=polyval(p,values);
LSB=abs(p(1))
s=abs(std(a2-f)/LSB)

[p,s]=polyfit(values,a3,1);
f=polyval(p,values);
LSB=abs(p(1))
s=abs(std(a3-f)/LSB)

figure(numer1);
plot(values,a2,'kv-',values,a3,'k*-',values,a1,'ko-');
grid on;
xlabel('DAC digital value');
ylabel('output current [uA]');
legend('R=220 kOhm, Vref=GND','R=3.7 MOhm, Vref=GND','R=3.7 MOhm, Vref=Vss',0);

figure(numer2);
plot(values,a2,values,a3,'kd-');
grid on;

cd /home/pawel/pliki/dane/chip1_plat;
channels=[0 1 2 3 4 5 8 12 57 58 59 60 61 62]

for i=1:length(channels)
  channel=num2str(channels(i));
  name=['plat',channel,'_chip1.dat']
  b0=importdata(name);
  b0=mean(b0');
  b(i,1:32)=b0;
  [p,s]=polyfit(values,b0,1);
  wzm(i)=abs(p(1));
  %f=polyval(p,values);
  %LSB=abs(p(1))
  %s=abs(std(a2-f)/LSB)
end













figure(100);
%plot(values,b1,values,b2);
