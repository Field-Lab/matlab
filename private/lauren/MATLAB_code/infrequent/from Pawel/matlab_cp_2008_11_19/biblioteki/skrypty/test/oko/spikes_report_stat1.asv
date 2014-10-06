nrchns=65;

s1=1000000;
dl=500000;
t=[1:dl];

samples=[s1 (s1+dl-1)];
channels=[6 13 20 29 36 65]; % +1 bo np. kanal 1 jest w pliku jako 2
%itd.
channel=29;
%channels=31;
name='Data009';
fs=20000;
read_param=struct('name',name,'header',206,'nrchns',nrchns,'channels',channel,'samples',samples);
s=readcnst(read_param);
time=[1:dl]/20;


signal=s(1:2:dl);

%sprawdzono, ze zmiana czestosci graqnicznej z 0.99 na np. 0.9 jest
%praktycznie bez znaczenia, dopiero przy 0.8 jest wyraznie gorzej
filtr_dl=20;
N=2;
[y,filtr]=oversampling(signal,N,filtr_dl,0.98);
y1=y(1,filtr_dl+1:filtr_dl+dl);

figure(1);
plot(t,s,t,y1-s)

figure(2);
for i=1:6
    read_param=struct('name',name,'header',206,'nrchns',nrchns,'channels',channels(1,i),'samples',samples);
    s=readcnst(read_param);
    signal=s(1:2:dl);
    [y,filtr]=oversampling(signal,N,filtr_dl,0.98);
    y1=y(1,filtr_dl+1:filtr_dl+dl);
    roznica=y1-s;
    figure(1);
    subplot(2,3,i);
    plot(t,s,t,roznica);
    figure(2);
    subplot(2,3,i);
    plot(roznica)
    std(roznica)
    figure(3);
    subplot(2,3,i);
    hist(roznica,100);
end