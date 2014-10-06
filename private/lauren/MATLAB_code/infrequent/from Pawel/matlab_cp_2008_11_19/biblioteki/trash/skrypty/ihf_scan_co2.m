cd /home/pawel/pliki/nauka/neuroplat;

ilf=20;
ihf=13;
igain=24;

a=importdata('pasmo_ch15_DAC_cal1_ilf20_scanihfco2_gain24_Vpol035_restfloating.dat');
wzmoc=zeros(16,35);
for i=1:16
	wzmoc(i,:)=a(i*65+19,:);
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

for i=1:16
	%if i~=10
		loglog(c1,wzmoc1(i,:));
		hold on;
	%else
	%	c2=[c1(1,1:24),c1(1,26:28)];
	%	s=[wzmoc1(i,1:24)),wzmoc1(1,26:28)];
	%	loglog(c2,s);
	%	hold on;
end
hold off;

grid on;
axis([10 10000 10 1000]);
%legend('external sinus','external square','internal calibration');
xlabel('frequency');
ylabel('gain');

figure(3)
for i=1:64
subplot(8,8,i);
s=a(650+i,:);
plot(s)
axis([0 35 0 0.3]);
end

