function [] = os1pluscz(fplus,wsp)
%fplus - częstotliwość dodawana do sygnału [Hz]


%odczyt
fidr=fopen('channel64.bin','r');
s1=fread(fidr,200000,'integer*2');
fclose(fidr);
dt=1/20000;
koniec=200000*dt;
t1=(0:dt:(koniec-dt))';
t2=1:200000;
s = s1 + 100*sin(2*pi*fplus*t1); %dodanie częstotliwości do sygnału

fs=fft(s);
figure(3)
subplot(2,1,1);
plot(abs(fs),'bd-')
%filtr notch
%fs = 20000;             %częstotliwość próbkowania
%f0 - częstotliwość notch
%f0=fplus;
%fn = fs/2;              %częstotliwość Nyquista
%freqRatio = f0/fn;      
%notchWidth = 0.1;       %szerokość notch

%Liczenie zer
%zeros = [exp( sqrt(-1)*pi*freqRatio ), exp( -sqrt(-1)*pi*freqRatio )];

%Liczenie biegunów
%poles = (1-notchWidth) * zeros;

%figure;
%zplane(zeros.', poles.');

%b = poly( zeros ); 
%a = poly( poles );

%figure(2);
%freqz(b,a,2000,fs)

%filtracja sygnału s
y = filter(wsp,s);
fy=fft(y);
subplot(2,1,2);
plot(abs(fy),'bd-')

x=s1-y;

figure (1);
subplot (3,1,1)
plot(t1,s)
xlabel('Sec')
axis tight
title('Oryginalny sygnał+dodana częstotliwość')
subplot (3,1,2)
plot(t1,y)
xlabel('Sec')

title('Sygnał przefiltrowany')
subplot(3,1,3)
plot(t1,s1)
title('Oryginalny sygnał')
xlabel('Sec')

end
