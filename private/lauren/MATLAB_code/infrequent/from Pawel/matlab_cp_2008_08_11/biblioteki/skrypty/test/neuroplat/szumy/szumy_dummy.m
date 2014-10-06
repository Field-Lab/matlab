% Skrypt do analizy wynikow z pomiarow szumow. Pomiar za pomoca progqramu w Labview "noise2a", lub "noise2", z katalogu "listopad2003/testy" w "moich dokumentach na Windows; UWAGA - plik zawiera dane przy pomiarze na kanale bocznym (ale zebrane za pomoca programu do czytania z multiplexera - bo mi sie nie chcialo przerabiac; zwykle probkowanie itd. - po prostu jako dane mozna wziac probki z dowolnego "kanalu"). 
%Plik zawiera wyestymowana gestosc widmowa. Fp=20000. N=32768. Normalizacja: FFT, modul kwadrat, razy 4, podzielic przez N, razy fp, razy 1e6. Aby obliczys sigme, trzeba wrocic do danych przed normalizacja, posumowac CALE widmo (a nie pol), spierwiastkowac, podzielic przez N.

cd /home/pawel/pliki/nauka/neuroplat;
a=importdata(' szumy_dummy_DAC_ilf20_ihf13_gain24_Vpol050.dat');
size(a)

fp=20000;
N=32768;
f=[1:16383]/N*fp;

figura1=1;
figura2=2;
figura3=3;

figure(figura1);
for i=1:64
	%figure(figura1);
	%subplot(8,8,i);
	%loglog(f,a(i+1,2:16384)*1e-6);
	%grid on;	
end

figure(figura3);
loglog(f,a(25,2:16384)*1e-6);
grid on;

figure(figura2);
plot(a(:,1))
grid on;
