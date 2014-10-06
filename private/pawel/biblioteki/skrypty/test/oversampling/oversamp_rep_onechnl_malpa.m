spikes=zeros(2,4);
thresholds=zeros(2,4);

N=2;
filtr_dl=20;
freq_gr=0.98;
filter_param=struct('N',N,'order',filtr_dl,'freq',freq_gr);

margins=[15 40];

figures=[10 11 12 13];


%cd I:\dane\2004\12maja\ForPawel\2003-08-06-0-016
cd /mnt/data4/dane/2004/12maja/ForPawel/2003-08-06-0-016

s1=importdata('electrode190.txt')';
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

%el_10: cb - 35,66 axonal - 64,197,207
%el167: cb - 1,25,141,281,242; axonal - 222,229,9,10,11,12,13,264
%el169: cb - 113,138; axonal - 1,2,285
%el190: cb - 244,290,295,20,80,162 (duuzo, b.wysokie); axonal - 126,129; dziwne -191,192,287,288
%el200: cb - 202;axonal - 8,26,89,92,224,244; dziwne: 54
%el353: cb - 196,197; axonal - 17,  dziwne:201,63,247,55
%el390: cb - 192,252,292 axonal - 105,186; dziwne - 26,68,151,248,51;
