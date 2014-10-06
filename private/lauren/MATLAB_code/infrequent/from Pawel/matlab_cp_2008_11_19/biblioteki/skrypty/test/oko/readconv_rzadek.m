clear;
cd /mnt/win3/oko/2000-12-12;

header=206;
nrchns=65;

start=150000;
dlugosc=1000;
koniec=start+dlugosc-1;
samples=[start koniec];

%matryca=[4 4];
%ilosc=;

dol=-800;
gora=800;

%start=2;
%channels=[start:(start+ilosc-1)];
channels=[6,21,23,34,41,55,65];
channels=[21:28];
channels=[21 65];
ilosc=length(channels);

figure(5);
r1=readconv('Data009conv',header,nrchns,channels,samples);
%8.8 MOhm
for i=1:(ilosc)
   subplot(ilosc,1,i);
   
   plot(r1(i,:));
   %axis([1 dlugosc dol gora]);
   %axis off;
end
