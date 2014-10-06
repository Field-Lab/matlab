cd /home/pawel/pliki/189e/pliki/oko/oko_64/marzec2002/
gain=importdata('chipA_all_gain.dat');
korr=importdata('chipA_all_korr.dat');
g=0.036;
tlumik=2.12/0.026;
f=[10 15 20 25 30 40 50 60 70 80 90 100 120 150 250 500 700 1000 1500 2000];

e=find(max(gain')>0.01);

wzm=gain/g*tlumik;
%for i=1:18
%wzm(:,i)=gain(:,i)/g;
%end

figure(122);
hold off;
clf(122);

for i=1:length(e)
%subplot(8,8,i);
loglog(f,wzm(e(i),:),'bd-')
%axis([10 2000 0 1400])
hold on;
grid on;
end
%grid on;

figure(123)
hold off;
clf(123);

for i=1:length(e)
subplot(8,8,i);
semilogx(f,korr(e(i),:),'bd-')
%axis([10 2000 0 1])
%hold on;
end
grid on;

'wdfwgwg'