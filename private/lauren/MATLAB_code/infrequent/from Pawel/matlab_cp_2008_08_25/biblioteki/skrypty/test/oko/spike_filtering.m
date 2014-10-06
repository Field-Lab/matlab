name='Data009';
fs=20000;

nrchns=65;
channels=65;
dl=4000;
s1=150000;
samples=[s1 (s1+dl-1)];
read_param=struct('name',name,'header',206,'nrchns',nrchns,'channels',channels,'samples',samples);
s=readcnst(read_param);


nrchns=65;
s1=3239700+134190;
dl=320;
samples=[s1 (s1+dl-1)];
channels=65;
%channels=31;
name='Data009';
fs=20000;
read_param=struct('name',name,'header',206,'nrchns',nrchns,'channels',channels,'samples',samples);
s1=readcnst(read_param);

s=s1;

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
   figure(1);
   subplot(2,3,i);
   plot(t,s,'bd-',t,y0,'g*-');
   axis([6.6 7.4 -1200 0]);
   xlabel('time [miliseconds]');
   legend('original signal',['cut-off freq. ' num2str(f(1,i)*10) 'kHz']);
   grid on;
end

figure(3);
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
   %axis([6 8 -1200 800]);
   legend('original signal', 'cut-off freq. 9 kHz','cut-off freq. 8 kHz','cut-off freq. 7 kHz','cut-off freq. 6 kHz','cut-off freq. 5 kHz','cut-off freq. 4 kHz')
end
grid on;
hold off;

