cd H:/pliki/nauka/stymulacja/chip/testy/2006_03_28/

rec=importdata('rec.dat');
art=importdata('art.dat');
artrec=importdata('art_and_rec.dat');

figure(1)
t=[1:4096];
%plot(t,rec(:,2))
plot(t,rec(:,2),t,art(:,2),t,artrec(:,2),t,artrec(:,2)-art(:,2));
grid on;
legend('rec','art','artrec','artrec-art');

a=artrec(:,2)-art(:,2);
a=a-mean(a);
b=rec(:,2);
b=b-mean(b);
figure(2)
plot(t,a,t,b);
grid on;

figure(3)
subplot(1,2,1);
plot(t,a);
grid on
subplot(1,2,2);
plot(t,b)
grid on

figure(4)
plot(t,a-b)
grid on

s=zeros(1:100,1);
for i=1:39
    s0=a(i*100+1:i*100+100,1);
    s=s+s0;
end
figure(6)
plot(s)