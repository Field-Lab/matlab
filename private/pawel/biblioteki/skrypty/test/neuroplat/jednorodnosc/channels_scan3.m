%Statystyka dla wszystkich kanalow - setup poprawiony. Pomiar dla wolnego zegara i w niskicz czestotliwosciach (patrz nazwa pliku).
%cd /home/pawel/pliki/nauka/neuroplat;
cd C:\pliki\nauka\neuroplat;

ilf=20;
ihf=13;
igain=24;
calib=1;

filename='pasmo_channelsscan_DAC_cal1_ilf01_ihf02_gain24_Vpol050_plytka1_AWG176MHz.dat';
%czestotliwosc AWG: 1.76 MHz, czyli jedna czwarta typowej
dzielnik_czest=4;
a1=importdata(filename);
wzmoc=zeros(64,35);

for i=1:64
	wzmoc(i,:)=a1(i*65+i+1,:);
end

wzmoc=wzmoc(:,2:29)*1e6/135;
c=importdata('czest.dat');
f=c(1,1:length(c))/dzielnik_czest;
clear c;

mnoznik=234/294; %z pomiaru:wysokosci piku oraz rzeczywistej amplitudy
wzmoc=wzmoc*mnoznik;

figure(1);
clf;

for i=1:64
	%if i~=8	& i~=61 & i~=1
		loglog(f,wzmoc(i,:));
		hold on;
	%end
end

grid on;
axis([40/dzielnik_czest 10000/dzielnik_czest 50 1000]);
xlabel('frequency');
ylabel('gain');
%print -depsc channels_scan


figure(2);
%wzm1=[wzmoc(1,18) wzmoc(2:64,15)'];
plot(max(wzmoc'),'bd-'); % bo......
grid on;
xlabel('channel number');
ylabel('gain');
text(60,950,['mean: ' num2str(mean(wzmoc(1:64,15)))]);
text(60,930,['mean: ' num2str(std(wzmoc(1:64,15))/mean(wzmoc(1:64,15))*100) '%']);

wzm1=wzmoc(1:64,15);

flow1=zeros(1,64);
fhigh1=zeros(1,64);

figure(4)
clf
for i=1:64
	max1=wzmoc(i,15);
	s=wzmoc(i,:);
	czest1=[20/dzielnik_czest:0.5:8000/dzielnik_czest];
	s=spline(f,s,czest1);
	a=find(s>max1/sqrt(2));
	flow1(1,i)=czest1(min(a));
	fhigh1(1,i)=czest1(max(a));
	%subplot(8,8,i);
	%if i~=8	& i~=61 & i~=1
		loglog(czest1,s);
	%end
	hold on;
end
hold off;
axis([40/dzielnik_czest 10000/dzielnik_czest 50 1000]);
grid on;

figure(3);
subplot(3,2,1);
hist_ph(wzm1);
xlabel('gain');
ylabel('nr of channels');

%figure(6)
subplot(3,2,2);
plot(wzm1,'bd-');
xlabel('channel');
ylabel('gain');
%axis([1 64 800 1000])
grid on;

%figure(4);
subplot(3,2,3);
hist_ph(flow1);
xlabel('low cut-off frequency');
ylabel('nr of channels');

subplot(3,2,4);
plot(flow1,'bd-');
xlabel('channel');
ylabel('low cut-off frequency');
grid on;
%axis([1 64 76 90]);
%axis([1 64])
%figure(5);

subplot(3,2,5);
hist_ph(fhigh1);
%axis([2400 2700 0 20])
xlabel('high cut-off frequency');
ylabel('nr of channels')

subplot(3,2,6);
plot(fhigh1,'bd-');
xlabel('channel');
ylabel('high cut-off frequency');
%axis([1 64 2400 2700]);
grid on;
