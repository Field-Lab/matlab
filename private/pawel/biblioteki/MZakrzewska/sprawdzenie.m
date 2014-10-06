function [] = sprawdzenie(fplus)
%fplus - cz�stotliwo�� dodawana do sygna�u [Hz]
%%%%odczyt
fidr=fopen('channel64.bin','r');
s1=fread(fidr,200000,'integer*2');
fclose(fidr);
dt=1/20000;
koniec=200000*dt;
t1=(0:dt:(koniec-dt))';
t2=1:200000;


%%%% stworzenie sygna�u zak��conego sinusem
mnoznik=100;
s=s1+mnoznik*sin(2*pi*fplus*t1);

%%%% Stworzenie filtru
fs=20000; %cz�stotliwo�� pr�bkowania  
f0=60/fs; %cz�stotliwo�� do usuni�cia znormalizowana
N=6; %liczba wsp�czynnik�w filtru
d  = fdesign.notch(N,f0,10,1);
    Hd = design(d);
    fvtool(Hd) %wizualizacja filtru (aby zmienic cz�stotliwo�c na znormalizowan� 
                %na wykresie odpowiedzi impulsowej nale�y wej�� w menu 
                %Analysis\Analysis Parameters i tam zaznaczy� normalised frequency)
wsp=Hd;

%%%% filtracja sygna�u s
tic
y=filter(wsp,s);
toc
%%%% Wykresy sygna�u
figure (1);
subplot(3,1,1)
plot(t1,s);
xlabel('Sec')
axis tight
title('Oryginalny sygna�+sinus 60')
subplot(3,1,2)
plot(t1,y,t1,s1);
xlabel('Sec')
title('Sygna� przefiltrowany')
h=gca;
set(h,'YLim', [-500 -300]);
subplot(3,1,3)
plot(t1,s1);
xlabel('Sec')
title('Oryginalny sygna�')
h=gca
set(h,'YLim', [-500 -300]);

%%%% Wykresy FFT
fs=fft(s); %fft sygna�u zak��conego
figure(2);
subplot(2,1,1);
plot(abs(fs),'bd-')
title('FFT sygna�u zak��conego')

fy=fft(y); %fft sygna�u przefiltrowanego
subplot(2,1,2);
plot(abs(fy),'bd-')



end
