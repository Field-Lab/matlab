cd H:\pliki\nauka\stymulacja\chip\testy\rozne;

a=importdata('sinus4_1uA_299kOhm.CSV');
t=a(:,1)*5/8*1000-28.8;  %skalowanie czasu - tak naprwde w pomiarze krok
%czasowy to 80us, tutaj skalowanie do 50us, zeby glupich pytan uniknac
s=(a(:,2)-mean(a(1:1000,2)))/299000*1e6;
a=plot(t,s,'k-');
set(a,'LineWidth',2);
axis([-0.5 5.5 -1.3 1.3]);
grid on;
h=gca;

%set(h,'XTickLabel',[])
%set(h,'YTickLabel',[])
fs=18;
set(h,'FontSize',fs);
set(h,'LineWidth',2)
a=xlabel('time [ms]');
set(a,'FontSize',fs);
a=ylabel('output current [ \muA]');
set(a,'FontSize',fs);