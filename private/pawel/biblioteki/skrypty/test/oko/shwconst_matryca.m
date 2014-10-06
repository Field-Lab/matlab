%lear;
clear;
%clf(2);
%cd /mnt/win3/oko/2000-12-11;
cd I:\dane\2000\2000-12-12;

header=206;
nrchns=65;
samples=[1 1600000];

matryca=[8 8];
ilosc=matryca(1)*matryca(2);

dol=-500;
gora=500;

start=2;
channels=[start:(start+ilosc-1)];
%channels=[2, 3, 5, 13, 14, 16, 21, 24:25,  27:28, 31:32, 34:36, 40, 44:53, 55, 62, 64];
ilosc=length(channels);

figure(5);
clf(5);
%r1=readconv('Data001conv',header,nrchns,channels,samples);
%8.8 MOhm
for i=1:(ilosc)
   r1=readconv('Data001conv',header,nrchns,channels(i),samples);
   b=shwconst2(r1,1000,2000,60);
   subplot(matryca(1),matryca(2),channels(i)-start+1);
   plot(b);
   axis([1 16000 -200 350]);
   %axis([samples(1) samples(2) dol gora]);
   %axis off;
end