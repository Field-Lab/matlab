nazwa='chipA_gr1_czest2.dat';
dane=importdata(nazwa);
d=extract_multi(dane,20,64);

figure(11);
clf(11);
for i=1:64
  subplot(8,8,i);
  plot(d(i,:))
  axis([0 32 -1.5 1]);
end
