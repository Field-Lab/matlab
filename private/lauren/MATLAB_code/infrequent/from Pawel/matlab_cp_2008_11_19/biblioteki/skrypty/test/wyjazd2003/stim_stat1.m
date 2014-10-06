%cd /home/pawel/pliki/dane/chip2_plat;
cd /home/pawel/pliki/189e/pliki/dane/chip2_plat;
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
  i 
  length(s1)
  s2=a2(i,length(s1):-1:1)/R*1000000;
  %t=[1:length(s1)];
  plot(t,s1,'k-',t,s2,'k--');
  if i==25 
    xlabel('Vstim [V]');
    ylabel('output current [uA]');
  end
  if i==1
    legend('increasing Vstim','decreasing Vstim',2);
  end
  text(-0.7,150,['DAC=' num2str(i-1)]);

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
  if i==25 
    xlabel('Vstim [V]');
    ylabel('output current [uA]');
  end

  if i==1
    legend('increasing Vstim','decreasing Vstim',2);
  end
  axis([-0.1 0.1 -20 20]);
  grid on;
  text(-0.07,15,['DAC=' num2str(i-1)]);
end

figure(205);
s1=a1(2,:)/R*1000000;
s2=a2(2,length(s1):-1:1)/R*1000000;
plot(t,s1,'k-',t,s2,'k--');
legend('increasing Vstim','decreasing Vstim',2);
grid on;
axis([-0.1 0.1 -2 2]);
h=gca;
set(h,'XTick',[-0.1 -0.05 0 0.05 0.1]);

cd /home/pawel/pliki/nauka/usa_2003/report/pictures;
obr_name='stim_stat';
print(203,'-depsc',obr_name);
print(203,'-dpng',obr_name);


obr_name='stim_stat_zoom';
print(204,'-depsc',obr_name);
print(204,'-dpng',obr_name);


obr_name='DAC_1_stim_stat';
print(205,'-depsc',obr_name);
print(205,'-dpng',obr_name);

