%Statystyka dla wszystkich kanalow - uwaga!! Kanal nr 1 daje bzdury dla typowego clocka, zmierzono go osobno z zegarem dwa razy wolniejszym. Stad te wszystki kwiatki przy rysowaniu, interpolowaniu splinem itd.

%cd /home/pawel/pliki/nauka/neuroplat;
cd C:\pliki\nauka\neuroplat;

ilf=20;
ihf=13;
igain=24;
calib=1;

a=importdata('pasmo_ch_scan_DAC_cal1_ilf20_ihf13_gain24_restfloating.dat');
%a=importdata('pasmo_scanChannels_DAC_cal1_ilf20_ihf13_gain24_Vpol035_restfloating.dat');
filename='pasmo_ch1_fp10kHz_DAC_cal1_ilf20_ihf13_gain24_restfloating.dat'
a1=importdata(filename);
wzmoc=zeros(64,35);

wzmoc(1,:)=a1(3,:); % bo osobny pomiar
for i=2:64
	wzmoc(i,:)=a(i*65+i+1,:);
end

wzmoc=wzmoc(:,2:29)*1e6/135;
c=importdata('czest.dat');
f=c(1,1:length(c));
clear c;
f1=f/2;

mnoznik=234/294; %z pomiaru:wysokosci piku oraz rzeczywistej amplitudy
wzmoc=wzmoc*mnoznik;

figure(1);
clf;

loglog(f(1,3:28)/2,wzmoc(1,3:28)); %bo osobny pomiar;
hold on;
for i=2:64
	if i~=8	& i~=61 & i~=1
		loglog(f,wzmoc(i,:));
		hold on;
	end
end

fontsize=18;
grid on;
axis([10 10000 50 1000]);
h=gca;
set(h,'FontSize',fontsize);
xlabel('frequency');
ylabel('gain');
print -depsc channels_scan



tekst1=['data path: wfitj71e/home/pawel/pliki/nauka/neuroplat'];
tekst2=['data file: ' filename];
tekst3='matlab path: wfitj71e/home/pawel/pliki/matlab/tulboksy/biblioteki/skrypty/test/neuroplat';
tekst4='matlab file: channels_scan1';
tekst5=['DACs: ilf='  num2str(ilf) ' ihf=' num2str(ihf) ' gain=scan' ' calib=' num2str(calib)];

tekstx=60;
teksty=300;
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
end %of "for"
end %of "if"



figure(2);
%wzm1=[wzmoc(1,18) wzmoc(2:64,15)'];
plot([wzmoc(1,18) wzmoc(2:64,15)'],'bd-'); % bo......
grid on;
xlabel('channel number');
ylabel('gain for 430 Hz');
text(60,950,['mean: ' num2str(mean(wzmoc(1:64,15)))]);
text(60,930,['mean: ' num2str(std(wzmoc(1:64,15))/mean(wzmoc(1:64,15))*100) '%']);

wzm1=wzmoc(1:64,15);

flow1=zeros(1,64);
fhigh1=zeros(1,64);

%dla kanalu 1 - ktory mierzony osobno (patrz importdata a1)
max1=wzmoc(1,18);
s=wzmoc(1,:);
czest1=[20:4000];
s=spline(f1,s,czest1);
a=find(s>max1/sqrt(2));
flow1(1,1)=czest1(min(a));
fhigh1(1,1)=czest1(max(a));

figure(4)
clf
for i=2:64
	max1=wzmoc(i,15);
	s=wzmoc(i,:);
	czest1=[20:0.5:4000];
	s=spline(f,s,czest1);
	a=find(s>max1/sqrt(2));
	flow1(1,i)=czest1(min(a));
	fhigh1(1,i)=czest1(max(a));
	%subplot(8,8,i);
	if i~=8	& i~=61 & i~=1
		loglog(czest1,s);
	end
	hold on;
end
axis([10 10000 50 1000]);

grid on;

figure(3);
subplot(3,2,1);
hist(wzm1);
axis([800 1000 0 20]);
xlabel('gain');
ylabel('nr of channels')
tx=930;
ty=13;
text(tx,ty,['mean:' num2str(mean(wzm1),4)]);
text(tx,ty-2,['sigma:' num2str(std(wzm1)/mean(wzm1)*100,'%3.1f') '%']);
grid on;

%figure(6)
subplot(3,2,2);
plot(wzm1,'bd-');
xlabel('channel');
ylabel('gain');
axis([1 64 800 1000])
grid on;

%figure(4);
subplot(3,2,3);
hist(flow1);
axis([76 90 0 20]);
xlabel('low cut-off frequency');
ylabel('nr of channels')
tx=84.5;
ty=18;
text(tx,ty,['mean:' num2str(mean(flow1),4)]);
text(tx,ty-2,['sigma:' num2str(std(flow1)/mean(flow1)*100,'%3.1f') '%']);
grid on;

subplot(3,2,4);
plot(flow1,'bd-');
xlabel('channel');
ylabel('low cut-off frequency');
grid on;
axis([1 64 76 90]);
%axis([1 64])
%figure(5);

subplot(3,2,5);
hist(fhigh1);
axis([2400 2700 0 20])
xlabel('high cut-off frequency');
ylabel('nr of channels')
tx=2580;
ty=18;
text(tx,ty,['mean:' num2str(mean(fhigh1),4)]);
text(tx,ty-2,['sigma:' num2str(std(fhigh1)/mean(fhigh1)*100,'%3.1f') '%']);
grid on;

subplot(3,2,6);
plot(fhigh1,'bd-');
xlabel('channel');
ylabel('high cut-off frequency');
axis([1 64 2400 2700]);
grid on;

%Offsety - czytane z innego pliku, zawierajacego zasadniczo gestosc widmowa szumow - pomiar odbywal sie ze zwartymi wejsciami, stad przy okazji zmierzyly sie offsety. Plik wykorzystywany w skrypcie szumy_Vpol_scan.m (Stan na 10 luty 2004).
a1=importdata('szumy_i_offsety_DAC_ilf20_ihf13_gain24_Vpol050.dat');

