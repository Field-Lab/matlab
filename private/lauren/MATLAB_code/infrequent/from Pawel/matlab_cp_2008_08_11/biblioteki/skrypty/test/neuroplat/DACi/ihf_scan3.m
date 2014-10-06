%Skrypt rysujacy rodzine charakterystyk neuroplata1.0.

cd /home/pawel/pliki/nauka/neuroplat;

figura1=1;
figura2=2;

filename='pasmo2_ch28_DAC_cal1_ilf02_ihfscan_gain24_Vpol050_chip1_AWG704MHz.dat';
%filename='pasmo2_ch28_DAC_cal1_ilf20_ihf13_gainscan_Vpol050_plytka1.dat';
a=importdata(filename);

ilf=20;
ihf=13;
gain=24;
AWGclock='7.04MHz';
AWGfile='CAL_SEQ.SEQ';
calib=1;
Vpol=0.5;


wzmoc=zeros(31,35);
for i=1:31
	wzmoc(i,:)=a(i*65+29,:);
	%wzmoc(i,:)=a(i*65+i+16,:);
	%wzmoc(i,:)=a(i*65+i+17,:);
end

wzmoc1=wzmoc(:,2:29)*1e6/135;
c=importdata('czest.dat');
c1=c(1,1:length(c));

mnoznik=234/294; %z pomiaru:wysokosci piku oraz rzeczywistej amplitudy
wzmoc1=wzmoc1*mnoznik;

figure(figura1);
clf;

for i=1:31
	%if i~=8 & i~=61
		loglog(c1,wzmoc1(i,:));
		%plot(c1,wzmoc1(i,:));
		hold on;
	%end
end
hold off;

fontsize=18;
grid on;
axis([10 10000 50 1000]);
h=gca;
set(h,'FontSize',fontsize);

%legend('external sinus','external square','internal calibration');
xlabel('frequency');
ylabel('gain');

tekstx=12;
teksty=200;

tekst1=['data path: wfitj71e/home/pawel/pliki/nauka/neuroplat'];
tekst2=['data file: ' filename];
tekst3='matlab path: wfitj71e/home/pawel/pliki/matlab/tulboksy/biblioteki/skrypty/test/neuroplat';
tekst4='matlab file: ihf_scan3';
tekst5=['DACs: ilf=' num2str(ilf) ' ihf=scan' ' gain=' num2str(gain) ' calib=' num2str(calib)];

opis=0;
if opis
for i=1:5
	switch i
	case 1
		tekst=tekst1;
	case 2
		tekst=tekst2;
	case 3
		tekst=tekst3;
	case 4
		tekst=tekst4;
	case 5
		tekst=tekst5;
	end
	ty=teksty*0.85^(i-1);
	a1=text(tekstx,ty,tekst);
	set(a1,'Interpreter','none');
	set(a1,'FontSize',12);
end
end

figure(figura2);
for i=1:64
	subplot(8,8,i);
	plot(a(65*27+i,:))
	axis([0 32 0 0.3]);
end



flow1=zeros(1,31);
fhigh1=zeros(1,31);

%dla kanalu 1 - ktory mierzony osobno (patrz importdata a1)
figure(3)
clf;
subplot(2,2,1);
for i=1:31
	max11=max(wzmoc1');
	max1=max11(i);
	s=wzmoc1(i,:);
	czest1=[12:0.5:8000];
	s=spline(c1,s,czest1);
	a=find(s>max1/sqrt(2));
	flow1(1,i)=czest1(min(a));
	fhigh1(1,i)=czest1(max(a));
	%subplot(8,8,i);
	%if i~=8	& i~=61 & i~=1
		loglog(czest1,s);
	%end
	hold on;
end
grid on;
axis([10 10000 50 1000]);
%axis([10 10000 50 1000]);
subplot(2,2,2);
plot(wzmoc1(:,18),'bd-');
axis([0 32 0 1000]);
grid on;
xlabel('ihf DAC');
ylabel('gain');
subplot(2,2,3);
plot(flow1,'bd-');
axis([0 32 12 20]);
grid on;
xlabel('ihf DAC');
ylabel('flow');
subplot(2,2,4);
plot(fhigh1,'bd-');
axis([0 32 0 5000]);
xlabel('ihf DAC');
ylabel('fhigh');
grid on;

figure(6);
plot(fhigh1,'bd-');
axis([0 32 0 5000]);
h=gca;
set(h,'FontSize',fontsize);
xlabel('ihf DAC [LSB]');
ylabel('fhigh [Hz]');
grid on;

tabelka=zeros(3,16);
tabelka(1,:)=max11(1:2:31);
tabelka(2,:)=flow1(1:2:31);
tabelka(3,:)=fhigh1(1:2:31);

fid=fopen('ihf_scan.txt','w');
fprintf(fid,'%6.1f %6.2f %6.0f\n',tabelka);
fclose(fid);

dane=[c1' wzmoc1'];
f=fopen('ihf_scan_data.txt','w');
fprintf(f,'%8.6f  %8.6f  %8.6f  %8.6f  %8.6f  %8.6f  %8.6f  %8.6f  %8.6f  %8.6f  %8.6f  %8.6f  %8.6f  %8.6f  %8.6f  %8.6f  %8.6f  %8.6f  %8.6f  %8.6f  %8.6f  %8.6f  %8.6f  %8.6f  %8.6f  %8.6f  %8.6f  %8.6f\n',dane);
%fwrite(f,dane,'double');
fclose(f);
q=importdata('ihf_scan_data.txt');
q
size(q)
size(dane)
