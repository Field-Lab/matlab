cd H:\pliki\nauka\stymulacja\chip\symulacje;
filename='15u_pos_96k.out';
bits=7;

a=importdata(filename);
x=a(:,1);
y=a(:,2);

[p,s]=polyfit(a(:,1),a(:,2),1);
w=polyval(p,a(:,1));
sum(w);
blad=(w-y)/max(abs(y))*2^bits;

[p2,Eteor2]=maxfit(x,y);
w4=polyval(p2,x);
blad4=(w4-y)/max(abs(y))*2^bits;
min(blad4)
max(blad4)
Eteor2*2^bits

figure(2);
plot(x,y,'b-',x,w,'k-',x,w4,'r-');
grid on;

figure(3);
plot(x,blad4,'bd-',x,blad,'kd-');
grid on;
