cd /home/pawel/pliki/nauka/neuroplat;

ilf=20;
ihf=13;
igain=24;

a=importdata('pasmo_ch_scan_DAC_cal1_ilf20_ihf13_gain24_restfloating.dat');
a=importdata('pasmo_scanChannels_DAC_cal1_ilf20_ihf13_gain24_Vpol035_restfloating.dat');
a=importdata('pasmo_ch15_DAC_cal1_ilf20_scanihf_gain24_Vpol035_restfloating.dat');
wzmoc=zeros(31,35);
for i=1:31
	wzmoc(i,:)=a(i*65+5,:);
	%wzmoc(i,:)=a(i*65+i+16,:);
	%wzmoc(i,:)=a(i*65+i+17,:);
end

wzmoc1=wzmoc(:,2:29)*1e6/135;
c=importdata('czest.dat');
c1=c(1,1:length(c));

mnoznik=234/294; %z pomiaru:wysokosci piku oraz rzeczywistej amplitudy
wzmoc1=wzmoc1*mnoznik;

figure(2);
clf;

for i=1:31
	%if i~=8 & i~=61
		loglog(c1,wzmoc1(i,:));
		hold on;
	%end
end
hold off;

grid on;
axis([10 10000 50 1000]);
%legend('external sinus','external square','internal calibration');
xlabel('frequency');
ylabel('gain');

