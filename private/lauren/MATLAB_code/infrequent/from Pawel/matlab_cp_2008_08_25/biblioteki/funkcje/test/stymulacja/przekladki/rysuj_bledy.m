function y=rysuj_bledy(current,res,bits);
%Funkcja rysuje przebiegi bledu fitowania. Aproksymacja dokonywana jest za
%pomoca funkcji maxfit. Funkcja wymaga istnienia w katalogu biezacym
%trzech plikow o nazwach typu: 1u_pos_20k.out, 1u_neg_20k.out,
%1u_input.out. 
%Wejscia:
%current - STRING, dla powyzszego przykladu '1u';
%res - STRING, dla powyzszego przykladu: '20k';
%bits - ilosc bitow, konieczna dla znormalizowania funkcji bledu i
%wyskalowania osi pionowej w LSB.

filename=[current '_neg_' res '.out'];
a=importdata(filename);
x=a(:,1);
y=a(:,2);
[p1,Eteor1]=maxfit(x,y);
w1=polyval(p1,x);
blad1=(w1-y)/max(abs(y))*2^bits;

filename=[current '_pos_' res '.out'];
a=importdata(filename);
x=a(:,1);
y=a(:,2);
[p2,Eteor2]=maxfit(x,y);
w2=polyval(p2,x);
blad2=(w2-y)/max(abs(y))*2^bits;

filename=[current '_input.out'];
a=importdata(filename);
x=a(:,1);
y=a(:,2);
[p3,Eteor3]=maxfit(x,y);
w3=polyval(p3,x);
blad3=(w3-y)/max(abs(y))*2^bits;

Eteor3*2^bits;

x=x/max(x);
plot(x,blad1,'b-',x,blad2,'k:',x,blad3,'r--');
grid on;

title(['range=' current ', res=' res]);
legend('neg. curr.','pos. curr.', 'input curr.');
xlabel('DAC output current [range]');
ylabel('linearity error [LSB]');