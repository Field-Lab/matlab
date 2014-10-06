function y=Rout1(current,res0,res1,rez,bits,osie);
%Funkcja rysuje stusunek pradu pos. do neg. jako funkcje wartosci DAC.
%Funkcja wymaga istnienia w katalogu biezacym
%dwoch plikow o nazwach typu: 1u_pos_20k.out, 1u_neg_20k.out,
%Wejscia:
%current - STRING, dla powyzszego przykladu '1u';
%res - STRING, dla powyzszego przykladu: '20k';
%bits - ilosc bitow, konieczna dla znormalizowania funkcji bledu i
%wyskalowania osi pionowej w LSB.

%1. czytanie danych z plikow
%a) dla pradu negatywnego (x - wartosci DACa, y - wartosci pradu)
%prad dla wyjscia zwartego:
filename=[current '_neg_' res0 '.out']
a=importdata(filename);
x_neg=a(:,1);
y_neg0=a(:,2);

%prad dla obciazenia:
filename=[current '_neg_' res1 '.out']
a=importdata(filename);
x_neg=a(:,1);
y_neg1=a(:,2);

%1b) dla pradu pozytywnego
%prad dla wyjscia zwartego:
filename=[current '_pos_' res0 '.out'];
a=importdata(filename);
x_pos=a(:,1);
y_pos0=a(:,2);

%prad dla obciazenia:
filename=[current '_pos_' res1 '.out'];
a=importdata(filename);
x_pos=a(:,1);
y_pos1=a(:,2);

Rout_neg=y_neg1./(y_neg0-y_neg1);
Rout_pos=y_pos1./(y_pos0-y_pos1);

y=Rout_pos;
value=x_pos/max(x_pos)*127;
a=semilogy(value,Rout_pos,'b-',value,Rout_neg,'r-');
set(a,'LineWidth',2);
grid on;
%xlabel('DAC current [LSB]');
%ylabel(['Rout [Rload=' res1 ']']);

%a=title([current 'A  R=' res1]);
a=text(70,1e4,[current 'A']);
set(a,'FontSize',16);
osie
%legend('neg. curr.','pos. curr.', 'input curr.');
if osie==1
    a=xlabel('DAC value [LSB]');
    set(a,'FontSize',16);
    a=ylabel('Rout [R]');
    set(a,'FontSize',16);
end
%y=a;
