%lear;
clear;
figura=12;
%clf(2);
%cd /mnt/win3/oko/2000-12-12;

%cd /home2/pawel/oko/2000-12-11; - linux
cd I:\dane\2000\2000-12-12; % windows

header=206;
%header=12;
nrchns=65;
samples=[200001 240000];
%st=8000000;
%samples=[st+1 st+200000-1];

matryca=[8 8];
ilosc=matryca(1)*matryca(2);
dlugosc=length(samples);

dol=-200;
gora=350;

start=2;
channels=[start:(start+ilosc-1)];
%channels=[22:31 33 54:56 60:65];
%channels=[23 26 56 61 63 65];
%channels=[2, 3, 5, 13, 14, 16, 21, 24:25,  27:28, 31:32, 34:36, 40, 44:53, 55, 62, 64];
ilosc=length(channels);

figure(figura);
clf(figura);
nazwa='Data009';
%read_param=2truct('name',nazwa,'header',206,'nrchns',65,'channels',channels(i),'samples',samples);
%r1=readcnst5('Data009conv',header,nrchns,channels,samples);
%8.8 MOhm
for i=1:(ilosc)
   read_param=struct('name',nazwa,'header',206,'nrchns',65,'channels',channels(i),'samples',samples);
   %r1=readcnst5(read_param);
   r1=readcnst(read_param);
   subplot(8,8,i);
   hand=gca;
   plot(r1);
   %axis([0 samples(2)-samples(1) -1000 1000]);
   %axis([samples(1) samples(2) dol gora]);
   %set(hand,'XTick',[0 6000 12000]);
   %set(hand,'XTickLabel',[0 5 10]);
   %set(hand,'YTickLabel',[]);
   %axis off;
end
