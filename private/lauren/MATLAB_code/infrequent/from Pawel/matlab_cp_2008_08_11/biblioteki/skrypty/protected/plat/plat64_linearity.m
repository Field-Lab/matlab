%cd /home/pawel/pliki/189e/pliki/dane/chip2_plat;
cd C:\pliki\189e\pliki\dane\chip2_plat;

R=10000;

a1=importdata('stim1_chip2.dat');
a2=importdata('stim2_chip2.dat');

numer1=203;
numer2=204;
t=[0:0.005:2]-1;

figure(numer1);
for i=1:32
  subplot(4,8,i);
  nr=num2str(i-1);
  tytul=['DAC value:' nr];
  title(tytul);
  s1=a1(i,:)/R*1000000;
  s2=a2(i,length(s1):-1:1)/R*1000000;
  %t=[1:length(s1)];
  plot(t,s1,'k-',t,s2,'k--');
  if i==1 
    xlabel('Vstim [V]');
    ylabel('output current [uA]');
  end
  if i==1
    legend('increasing Vstim','decreasing Vstim',2);
  end
  
  axis([-1 1 -200 200]);
  grid on;
end

figure(numer2);
for i=1:32
  subplot(4,8,i);
  s1=a1(i,:)/R*1000000;
  s2=a2(i,length(s1):-1:1)/R*1000000;
  %t=[1:length(s1)];
  plot(t,s1,'k-',t,s2,'k--');
  if i==1
    legend('increasing Vstim','decreasing Vstim',2);
  end
  axis([-0.1 0.1 -0.8 0.8]);
  grid on;
end
