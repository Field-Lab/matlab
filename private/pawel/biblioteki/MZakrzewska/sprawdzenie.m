function [] = sprawdzenie(fplus)
%fplus - czêstotliwoœæ dodawana do sygna³u [Hz]
%%%%odczyt
fidr=fopen('channel64.bin','r');
s1=fread(fidr,200000,'integer*2');
fclose(fidr);
dt=1/20000;
koniec=200000*dt;
t1=(0:dt:(koniec-dt))';
t2=1:200000;


%%%% stworzenie sygna³u zak³óconego sinusem
mnoznik=100;
s=s1+mnoznik*sin(2*pi*fplus*t1);

%%%% Stworzenie filtru
fs=20000; %czêstotliwoœæ próbkowania  
f0=60/fs; %czêstotliwoœæ do usuniêcia znormalizowana
N=6; %liczba wspó³czynników filtru
d  = fdesign.notch(N,f0,10,1);
    Hd = design(d);
    fvtool(Hd) %wizualizacja filtru (aby zmienic czêstotliwoœc na znormalizowan¹ 
                %na wykresie odpowiedzi impulsowej nale¿y wejœæ w menu 
                %Analysis\Analysis Parameters i tam zaznaczyæ normalised frequency)
wsp=Hd;

%%%% filtracja sygna³u s
tic
y=filter(wsp,s);
toc
%%%% Wykresy sygna³u
figure (1);
subplot(3,1,1)
plot(t1,s);
xlabel('Sec')
axis tight
title('Oryginalny sygna³+sinus 60')
subplot(3,1,2)
plot(t1,y,t1,s1);
xlabel('Sec')
title('Sygna³ przefiltrowany')
h=gca;
set(h,'YLim', [-500 -300]);
subplot(3,1,3)
plot(t1,s1);
xlabel('Sec')
title('Oryginalny sygna³')
h=gca
set(h,'YLim', [-500 -300]);

%%%% Wykresy FFT
fs=fft(s); %fft sygna³u zak³óconego
figure(2);
subplot(2,1,1);
plot(abs(fs),'bd-')
title('FFT sygna³u zak³óconego')

fy=fft(y); %fft sygna³u przefiltrowanego
subplot(2,1,2);
plot(abs(fy),'bd-')



end
