cd pliki/oko/oko_64/marzec2002/

gain=importdata('chipA_all_gain_nast2.dat');
korr=importdata('chipA_all_korr_nast2.dat');
g=0.036;
tlumik=2.12/0.026;
f=[10 15 20 25 30 40 50 60 70 80 90 100 120 150 250 500 700 1000 1500 2000];

e=find(max(gain')>0.3);

'sgggggggggggggggg'

wzm=gain/g*tlumik;
%for i=1:18
%wzm(:,i)=gain(:,i)/g;
%end

figure(22);
hold off;
clf(22);

for i=1:32
%subplot(2,1,1);
%semilogx(f,wzm(e(2*i-1),:),'bd-')
loglog(f,wzm(e(2*i-1),:),'bd-')
axis([10 2000 0 1400])
hold on;
%subplot(2,1,2);
%semilogx(f,wzm(e(2*i),:),'g*-')
loglog(f,wzm(e(2*i),:),'g*-')
axis([10 2000 0 1400])
hold on;
end
%grid on;

figure(23)
hold off;
clf(23);

for i=1:length(e)
subplot(8,8,i);
semilogx(f,korr(e(i),:),'bd-')
axis([10 2000 0 1])
%hold on;
end
grid on;
