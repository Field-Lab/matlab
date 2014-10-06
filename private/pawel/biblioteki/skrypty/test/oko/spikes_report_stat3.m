%cd /mnt/data4/dane/2000/2000-12-12; %linux
cd I:\dane\2000\2000-12-12 %windows

nrchns=65;

s1=1000000;
s1=s1+1500000;
dl=1000000;
t=[1:dl];

samples=[s1 (s1+dl-1)];
channels=[6 13 20 29 36 65]; % +1 bo np. kanal 1 jest w pliku jako 2
%itd.
channel=29;
%channels=31;
name='Data009';
fs=20000;

%sprawdzono, ze zmiana czestosci graqnicznej z 0.99 na np. 0.9 jest
%praktycznie bez znaczenia, dopiero przy 0.8 jest wyraznie gorzej
filtr_dl=20;
N=2;

marg_left=10; %tyle probek "na lewo" od poczatku spike'a zalicza sie do spike'a
marg_right=50;

detekcja=[400,100,100,200,100,300];

figure(3);
hold off;
clf;
figure(4);
hold off;
clf;
figure(11);
hold off;
clf;
figure(5);
hold off;
clf;
figure(6);
hold off;
clf;


ilosc=zeros(1,6);

kolory={'bd','gd','rd','cd','md','kd'};

channels=[2 3 5 6 7 9 12 13 14 16 18:25 27:33 35 37 38 39 41 50 51:54 58 59 61 63 64];
l_channels=length(channels);

clear en blad wynik wynik0;
%for i=1:l_channels
%    read_param=struct('name',name,'header',206,'nrchns',nrchns,'channels',channels(i)+1,'samples',samples);
%    s=readcnst(read_param);
%    figura=5;
%    figure(figura);
%   subplot(8,5,i);
%    y=blad_vs_energy(s,figura);
    %hold on;
    %end
%hold off;

channels=


figures=[1 2 3 4];
for i=1:1
    read_param=struct('name',name,'header',206,'nrchns',nrchns,'channels',channels(i)+1,'samples',samples);
    s=readcnst(read_param);    
    figura=6;
    figure(figures(i));
    hold off;
    clf;
    y=oversamp_details1(s,figura,marg_left,marg_right);
end
hold off;