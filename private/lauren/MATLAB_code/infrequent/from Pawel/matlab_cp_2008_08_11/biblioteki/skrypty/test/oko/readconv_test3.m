clear;
cd /mnt/win3/oko/2000-12-12;

header=206;
nrchns=65;
samples=[1 20000];

matryca=[8 8];
ilosc=matryca(1)*matryca(2);

dol=0;
gora=500;

start=2;
%channels=[start:(start+ilosc-1)];
channels=[2:8, 13:16, 20:25, 27:32, 34:37, 40, 44:55, 62, 64];
ilosc=length(channels)

figure(1);
r1=readconv('Data008conv',header,nrchns,channels,samples);
%8.8 MOhm
for i=1:ilosc
   %subplot(matryca(1),matryca(2),channels(i));
   a=r1(i,:);
   srednia(i)=mean(a);
   poziom(i)=std(a);
   %plot(r1(i,:));
   %axis([samples(1) samples(2) dol gora]);
end
plot(srednia);
figure(2);
plot(poziom);


figure(3);
%30 MOhm
r2=readconv('Data010conv',header,nrchns,channels,samples);
for i=1:ilosc
   subplot(matryca(1),matryca(2),channels(i));
   plot(r2(i,:));
   axis([samples(1) samples(2) dol gora]);   
end

figure(5);
%130 MOhm
r3=readconv('Data012conv',header,nrchns,channels,samples);
for i=1:ilosc
   subplot(matryca(1),matryca(2),channels(i));
   plot(r3(i,:));
   axis([samples(1) samples(2) dol gora]);   
end

ch=[2:8, 13:16, 20:25, 27:32, 34:37, 40:41, 44:55, 62, 64];