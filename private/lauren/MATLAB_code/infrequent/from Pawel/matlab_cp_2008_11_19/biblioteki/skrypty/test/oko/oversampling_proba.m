readcnst_channel_malezmiany;

s=s65;
signal=s(1:2:dl);
dl=length(s);

filtr_dl=16;
N=2;
[y,filtr]=oversampling(signal,N,filtr_dl,0.99);
y=N*y;
y0=y(1,1+filtr_dl:length(s)+filtr_dl);

%y=N*y;
figure(3);
subplot(2,2,1);
plot(s,'bd-');
grid on;
axis([175 200 -1500 1000]);
subplot(2,2,2);
plot(abs(fft(s)));
axis([0 round(length(s)/2) 0 5000]);
subplot(2,2,3);
plot(y0,'bd-');
grid on;
axis([175 200 -1500 1000]);
subplot(2,2,4);
plot(abs(fft(y0)));
axis([0 round(length(s)/2) 0 5000]);

figure(9);
t=[1:length(s)];
plot(t,s,'bd-',t,y0,'gd-');

size(filtr)

smax=max(abs(signal))
ymax=max(abs(y))

find(abs(s)>smax-1)
find(abs(y)>ymax-1)

%y0=y(1,1+filtr_dl:length(s)+filtr_dl);
figure(5);
subplot(2,2,1);
blad=y0-s;
plot(blad);
subplot(2,2,2);
plot(abs(fft(blad)));
blad0=blad(1,1:500);
subplot(2,2,3);
plot(blad0,'bd-');
subplot(2,2,4);
plot(abs(fft(blad0)));

std(blad)

spike1=s(1,1:500);
figure(6);
subplot(2,1,1);
plot(spike1);
subplot(2,1,2);
plot(abs(fft(spike1)));

figure(7)
subplot(2,1,1);
plot(filtr);
subplot(2,1,2);
plot(abs(fft(filtr,1000)));

