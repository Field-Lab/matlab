cd /mnt/win3/oko/2000-12-11/;

%channels=[9 10 11 22 27 52:56 60:65];
%channels=[27 52 53 55 61:65];
%channels=[58:60];
%channels=[22:31 33 54:56 60:65];
%channels=[2:65];
channels=[32 33];
name='Data002conv';

dlugosc=400000;
start=[1 2500000 6000000 14000000];

N=256;
M=250;
fp=20000;
lsb=1;

figura1=1;
figura2=4;
figura3=5;

%hannels=[60 61];
%channels=[26 56 63 65];

w=shwnoise(name,channels,dlugosc,start,N,M,fp,lsb,1,1,5);

name='Data001conv';
%w=shwnoise(name,channels,dlugosc,start,N,M,fp,lsb,1,1,3);

name='Data002conv';
%w=shwnoise(name,channels,dlugosc,start,N,M,fp,lsb,1,1,4);
