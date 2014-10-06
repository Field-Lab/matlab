cd /home/pawel/pliki/nauka/neuroplat;

ilf=20;
ihf=13;
igain=24;

a=importdata('pasmo_ch_scan_DAC_cal1_ilf20_ihf13_gain24_restfloating.dat');
a=importdata('pasmo_scanChannels_DAC_cal1_ilf20_ihf13_gain24_Vpol035_restfloating.dat');
wzmoc=zeros(64,35);
for i=1:64
	wzmoc(i,:)=a(i*65+i+1,:);
end

wzmoc1=wzmoc(:,2:29)*1e6/135;
c=importdata('czest.dat');
c1=c(1,1:length(c));

mnoznik=234/294; %z pomiaru:wysokosci piku oraz rzeczywistej amplitudy
wzmoc1=wzmoc1*mnoznik;

figure(1);
clf;

for i=1:64
	if i~=8	& i~=61
		loglog(c1,wzmoc1(i,:));
		hold on;
	end
end
hold off;

grid on;
axis([10 10000 50 1000]);
%legend('external sinus','external square','internal calibration');
xlabel('frequency');
ylabel('gain');

figure(3);
plot(wzmoc1(:,15),'bd-');
grid on;
xlabel('channel number');
ylabel('gain for 430 Hz');
text(60,950,['mean: ' num2str(mean(wzmoc1(1:64,15)))]);
text(60,930,['mean: ' num2str(std(wzmoc1(1:64,15))/mean(wzmoc1(1:64,15))*100) '%']);

figure(4);
for i=1:28
	subplot(5,6,i);
	plot(wzmoc1(:,i),'bd-');
	grid on;
	axis([0 65 80 800]);
end

