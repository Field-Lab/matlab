spikes=zeros(2,4);
thresholds=zeros(2,4);

N=2;
filtr_dl=20;
freq_gr=0.98;
filter_param=struct('N',N,'order',filtr_dl,'freq',freq_gr);

margins=[15 40];

figures=[10 11 12 13];

cd I:\dane\2004\26maja\2003-12-08-0-001;
%cd I:\dane\2004\26maja\2004-01-13-0-001;
%cd /mnt/data4/dane/2004/26maja/2003-12-08-0-001;
s1=importdata('electrode45.txt')';
s1=s1-mean(s1);

prog=60;
histereza=30;
znak=-1;
detect_param=struct('prog',prog,'histereza',histereza,'znak',znak);
%figure(figura1);

[y1,y2]=inveyefilter_hayes(200,50,50,3100,20000); %oszacowanie efektu filtrowanai w neurochipie
figure(156);
y1bis=[y1(1,202:length(y1)) y1(1,1:201)]; 
f_y1=fft(y1bis);
freq_y1=[0:length(y1)-1]/length(y1)*20000;
subplot(4,2,1);
semilogx(freq_y1,abs(f_y1));
h=gca;
set(h,'XLim',[50 10000]);
grid on;
subplot(4,2,3);
semilogx(freq_y1,angle(f_y1),'bd-');
h=gca;
set(h,'XLim',[50 10000]);
grid on;
subplot(4,2,5);
plot(y1);
%grid on;

f_y2=fft(y2);
freq_y1=[0:length(y1)-1]/length(y1)*20000;
subplot(4,2,2);
semilogx(freq_y1,abs(f_y2));
h=gca;
set(h,'XLim',[50 10000]);
grid on;
subplot(4,2,4);
semilogx(freq_y1,angle(f_y2),'bd-');
h=gca;
set(h,'XLim',[50 10000]);
grid on;
subplot(4,2,6);
plot(y2);
%grid on;

subplot(4,2,7);
plot(conv(y1,y2));


%y2 to filtr odwrotny, jest on potem przekazywany do funkcji robioacej
%statystyki.
figure(150);
s=s1(1,1:50000);
s=s-mean(s);
sy=conv(s,y2);
size(sy)
sy=sy(201:length(s)+200);
sy=sy-mean(sy);
t=[1:length(s)]/20;
subplot(2,1,1);
plot(s);
axis([6400 7000 -800 200]);
grid on;
subplot(2,1,2);
plot(sy);
axis([6400 7000 -800 200]);
grid on;
figure(253);
plot(t,s,'b-',t,sy,'r-');
grid on;

%break;

y=spikes_oversamp_stat1(s1,detect_param,filter_param,margins,figures,y2);
%I:\dane\2004\26maja\2003-12-08-0-001 - kolejne pliki:
%electrode45: cell body: 85, axonal: 185,300
%electrode137: cell body: 30, 240 axonal:231, 234
%electrode148: cell body 227, axonal: 106, 254
%electrode155: 2,63 axonal: 285
%electrode169: cb: 158, 181 axonal: 5,285
%el221: cb - 10,147 axonal:none
%el298: cb - 1,20; axonal:227, 228
%el362: cb - 4,5; axonal: 3, 225, 227, 281 (mnooostwo!)
%el416: cb - 61,68,167; axonal: 80, 161, 180 (duzo)
%el499: cb - 14, 17; axonal:none

%I:\dane\2004\26maja\2004-01-13-0-001
%electrode35: nic tam nie ma
%electrode130: cb: 4,9,11; axonal: 281, 294 (duuuuzo)
%el152: 181, 182, 183 (mnoostwo), axonal:none
%el313: nic ciekawego
%el337: cb - 2,3 (b.niskie), axonal:none
%el370: nic
%el424: nic
%el471: nic
%el474: nic
