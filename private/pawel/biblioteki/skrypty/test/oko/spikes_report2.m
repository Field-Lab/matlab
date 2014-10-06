%Skrypt dokonujacy:
%-plotowania okreslonego spikea i jego widma;
%-plotowania widm grupy spikow;
%-filtrowania wybranego fragmentu sugnali i jego prezeentacji

cd /home2/pawel/oko/2000-12-12;

s1=150000;
dl=4000;

nrchns=65;
%s1=3239700+134190;
%dl=320;
samples=[s1 (s1+dl-1)];
channels=65;
%channels=31;
name='Data009';
fs=20000;
read_param=struct('name',name,'header',206,'nrchns',nrchns,'channels',channels,'samples',samples);
s=readcnst(read_param);
time=[1:dl]/20;


%nrchns=65;
%channels=65;
%dl=500;
%s1=149950;
%samples=[s1 (s1+dl-1)];
%read_param=struct('name',name,'header',206,'nrchns',nrchns,'channels',channels,'samples',samples);
%s1=readcnst(read_param);
%time=[1:dl]/20000*1000;

%Seria spikeow - sygnal:
figure(1);
subplot(2,2,1);
plot(time,s);
axis([0 200 -1200 800]);
xlabel('time [miliseconds]');
ylabel('signal value');
text(162,760,'2000/12/12');
text(162,695,'Data009');
text(162,630,'channel 64');
grid on;
title('Fig. 1: spikes burst');

%Jeden wyciety spike - sygnal:
figure(1);
subplot(2,2,3);
plot(time,s);
axis([8.3 12 -1200 700]);
axis([6 14 -1200 800]);
xlabel('time [miliseconds]');
ylabel('signal value');
grid on;
title('Fig. 3: spike number one from figure 1 (look on the time base)');

%Tenze spike - widmo:
figure(1);
subplot(2,2,4);
s0=s(1,120:280);
dl0=length(s0);
f=fft(s0)/dl0*2*dl0/90;
cz=[0:dl0-1]/dl0*fs;
plot(cz,abs(f));
axis([0 fs/2 0 300]);
xlabel('frequency [Hz]');
ylabel('abs. value of Fourier coefficient');
grid on;
title('Fig. 4: Fourier Transform of the spike from figure 3');

f0=f;
cz0=cz;

%Widmo serii z rys. 1
figure(1);
subplot(2,2,2);
f=fft(s)/dl*2;
cz=[0:dl-1]/dl*fs;
plot(cz,abs(f));
axis([0 fs/2 0 50]);
xlabel('frequency [Hz]');
ylabel('abs. value of Fourier coefficient');
grid on;
title('Fig. 2: Fourier Transfrom of the spike burst (signal from figure 1)');


dl=length(s);
t=[1:dl]/fs*1000;

filtr_dl=16;

f=[0.4 0.5 0.6 0.7 0.8 0.9];
for i=1:6
   filtr=fir1_ph(filtr_dl,f(1,i));
   size(filtr)
   y=conv(s,filtr);
   size(y)
   y0=y(1,1+filtr_dl/2:dl+filtr_dl/2);
   %figure(3);
   %subplot(2,3,i);
   %plot(y);i
   figure(5);
   subplot(2,3,i);
   plot(t,s,'bd-',t,y0,'g*-');
   axis([8.9 9.5 -1200 0]);
   xlabel('time [miliseconds]');
   legend('original signal',['cut-off freq. ' num2str(f(1,i)*10) 'kHz']);
   grid on;
end
subplot(2,3,1);
title('Fig 6: falling slope of the filtered spike from figure 3');

figure(6);
clf
plot(t,s);
hold on;
styles=['b-' 'g-'];
for i=1:6
   filtr=fir1_ph(filtr_dl,f(1,7-i));
   size(filtr)
   y=conv(s,filtr);
   size(y)
   y0=y(1,1+filtr_dl/2:dl+filtr_dl/2);
   %subplot(6,2,2*i);
   plot(t,y0);
   axis([6 14 -1200 800]);
   legend('original signal', 'cut-off freq. 9 kHz','cut-off freq. 8 kHz','cut-off freq. 7 kHz','cut-off freq. 6 kHz','cut-off freq. 5 kHz','cut-off freq. 4 kHz')
   xlabel('time [miliseconds]');
   ylabel('signal value');
end
grid on;
title('Fig. 5: filtered spike from figure 3 - for various cut-off frequencies');
hold off;


%Nadprobkowanie:
signal=s(1:2:dl);

filtr_dl=36;
N=2;
[y,filtr]=oversampling(signal,N,filtr_dl,0.99);
%y=N*y;
y0=y(1,1+filtr_dl:length(s)+filtr_dl);
figure(7);
subplot(2,2,1);
plot(t,s,'k-',t,y0,'k--');
axis([6 14 -1200 800]);
legend('original signal','interpolated signal');
xlabel('time [miliseconds]');
ylabel('signal value');
grid on;
title('Fig. 7: original signal and simulation of 10 kHz ADC sampling with software Nyquist interpolation');

figure(7);
subplot(2,2,2);
plot(t,s,'bd-',t,y0,'g*-');
axis([8.9 9.5 -1200 0]);
legend('original signal','interpolated signal');
xlabel('time [miliseconds]');
ylabel('signal value');
grid on;
title('Fig. 8: the zoom of the signals from Fig. 7'); 

figure(7);
subplot(2,2,3);
y1=y0(1,120:280);
f1=fft(y1)*2/90;
plot(cz0,abs(f0),'k-',cz0,abs(f1),'k--');
axis([0 10000 0 300]);
xlabel('frequency [Hz]');
ylabel('abs. value of Fourier coefficient');
legend('original signal','interpolated signal');
title('Fig.9: Fourier Transforms of the signals from Fig. 7');
grid on;

figure(7);
subplot(2,2,4);
df=20;
czest=[0:20:19980];
plot(czest,abs(fft(filtr,1000)));
axis([0 10000 0 1.2]);
xlabel('frequency [Hz]');
ylabel('abs. value of Fourier coefficient');
title('Fig 10: Fourier Transform of the used Nyquist interpolation filter');
grid on;

