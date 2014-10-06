%lear;
clear;
figura=3;
%clf(2);
cd /mnt/win3/oko/2000-12-11;

header=206;
nrchns=65;
samples=[1 50000];

matryca=[8 8];
ilosc=matryca(1)*matryca(2);

dol=-500;
gora=500;

start=2;
channels=[start:(start+ilosc-1)];
%channels=[2, 3, 5, 13, 14, 16, 21, 24:25,  27:28, 31:32, 34:36, 40, 44:53, 55, 62, 64];
ilosc=length(channels)

figure(figura);
clf(figura);
r1=readconv('Data001conv',header,nrchns,channels,samples);
%8.8 MOhm
for i=1:(ilosc)
   subplot(matryca(1),matryca(2),channels(i)-start+1);
   plot(hist(r1(i,:)));
   %axis([samples(1) samples(2) dol gora]);
   %axis off;
end
