clear;

header=206;
nrchns=65;
samples=[1 20000];

%channels=[start:(start+ilosc-1)];
channels=[2, 3, 5, 13, 14, 16, 21, 24:25,  27:28, 31:32, 34:36, 40, 44:53, 55, 62, 64];
ilosc=length(channels)

C=150*1e-12;
R1=8.8*1e+6;
R2=30*1e+6;
R3=130*1e+6;

f1=13;
f2=3000;

rzad_filtru=5000;

%figure(4);
r1=readconv('Data008conv',header,nrchns,channels,samples);  %8.8 MOhm
%f3=1/(C*R1);
%r11=inveyefilter(rzad_filtru,f1,f2,f3,20000);

r2=readconv('Data010conv',header,nrchns,channels,samples);  %30 MOhm
r3=readconv('Data012conv',header,nrchns,channels,samples);	%130 MOhm

[sr1,od1]=channel_stat(r1);
[sr2,od2]=channel_stat(r2);
[sr3,od3]=channel_stat(r3);

a=size(r1);
b=[1:a(1)];

figure(4);
plot(b,sr1,b,sr2,b,sr3);  %niebieska, zielona, czerwona
figure(5);
plot(b,od1,'bd-',b,od2,'gd-',b,od3,'rd-');

%8.8 MOhm
%ch=[2:8, 13:16, 20:25, 27:32, 34:37, 40:41, 44:55, 62, 64];