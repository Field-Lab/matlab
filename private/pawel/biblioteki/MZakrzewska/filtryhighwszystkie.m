tmin=3.8;
tmax=5.2;
vmin=-0.5;
vmax=0.3;

figure(1)
clf
plot (t1, s1)
axis([tmin tmax  vmin vmax])
xlabel ('t[sec]')
ylabel('amplitude')

fs=10000;
F=0:1/fs:1;
H2=[zeros(1,20+1),ones(1,fs-20)];
h2=fir2(200,F,H2);
%figure
%freqz(h2,1);
y200=fftfilt(h2,s1);
figure(2)
clf
plot ((t1+0.004),s1,'b','LineWidth',1.5)
hold on
plot (t1,y200,'r','LineWidth',1.5)
axis([tmin tmax  vmin vmax])
xlabel ('t[sec]')
ylabel('amplitude')

F=0:1/fs:1;
H2=[zeros(1,35+1),ones(1,fs-35)];
h2=fir2(200,F,H2);
%figure
%freqz(h2,1);
y300=fftfilt(h2,s1);
figure(3)
clf
plot ((t1+0.004),s1,'b','LineWidth',1.5)
hold on
plot (t1,y300,'r','LineWidth',1.5)
axis([tmin tmax  vmin vmax])
xlabel ('t[sec]')
ylabel('amplitude')


F=0:1/fs:1;
H2=[zeros(1,60+1),ones(1,fs-60)];
h2=fir2(200,F,H2);
%figure
%freqz(h2,1);
y600=fftfilt(h2,s1);
figure(4)
clf
plot ((t1+0.004),s1,'b','LineWidth',1.5)
hold on
plot (t1,y600,'r','LineWidth',1.5)
axis([tmin tmax  vmin vmax])
xlabel ('t[sec]')
ylabel('amplitude')

F=0:1/fs:1;
H2=[zeros(1,100+1),ones(1,fs-100)];
h2=fir2(200,F,H2);
%figure
%freqz(h2,1);
y1000=fftfilt(h2,s1);
figure(5)
clf
plot ((t1+0.004),s1,'b','LineWidth',1.5)
hold on
plot (t1,y1000, 'r','LineWidth',1.5)
axis([tmin tmax  vmin vmax])
xlabel ('t[sec]')
ylabel('amplitude')


figure(6)
clf
plot (t1,y200,'r','LineWidth',3)
hold on
plot(t1,y300,'g','LineWidth',3)
hold on
plot(t1, y600,'y','LineWidth',3)
hold on
plot(t1, y1000,'c','LineWidth',3)
hold on
plot ((t1+0.004),s1,'b','LineWidth',3)
axis([tmin tmax  vmin vmax])
xlabel ('t[sec]')
ylabel('amplitude')
