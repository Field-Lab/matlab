%wczytanie danych
fidr=fopen('channel64.bin','r');
s1=fread(fidr,200000,'integer*2');
fclose(fidr);

%stworzenie wektora czasu
dt=1/20000;
koniec=200000*dt;
t1=(0:dt:(koniec-dt))'; %wektor czasu [s]
t2=1:200000;


s=s1;%+sin(2*pi*60)  s - sygna�


fs = 20000;             % cz�stotliwo�� pr�bkowania
f0=200;                 % cz�stotliwo�� graniczna


F=0:1/10000:1;
H2=[zeros(1,f0+1),ones(1,10000-f0)];
h2=fir2(200,F,H2);
%figure
%freqz(h2,1)
y1=fftfilt(h2,s);


fcut=60;               % cz�stotliwo�� notch
fn = fs/2;              %cz�stotliwo�� Nyquista
freqRatio = fcut/fn;      
notchWidth = 0.1;       %szeroko�� notch

%Liczenie zer
zeros = [exp( sqrt(-1)*pi*freqRatio ), exp( -sqrt(-1)*pi*freqRatio )];

%Liczenie biegun�w
poles = (1-notchWidth) * zeros;

b = poly( zeros ); 
a = poly( poles );

%filtracja sygna�u x
y2 = filter(b,a,y1);
figure
% pierwszy wykres
subplot(2,1,1)
plot(t1,s1);
xlabel('Sec')
axis tight
title('Oryginalny sygna�')
%drugi wykres
subplot(2,1,2)
plot(t1,y2);
xlabel('kHz')
title('Sygna� przefiltrowany')