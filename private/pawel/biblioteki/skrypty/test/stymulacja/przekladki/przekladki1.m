cd H:\pliki\nauka\stymulacja\chip\symulacje;
filename='60u_pos.out';
bits=7;

a=importdata(filename);
x=a(:,1);
y=a(:,2);
[p,s]=polyfit(a(:,1),a(:,2),1);
p
w=polyval(p,a(:,1));
sum(w);
blad=(y-w)/max(abs(y))*2^bits;

figure(1);
plot(blad);
grid on;

figure(2);
plot(x,y,'bd-',x,w,'kd-');
grid on;