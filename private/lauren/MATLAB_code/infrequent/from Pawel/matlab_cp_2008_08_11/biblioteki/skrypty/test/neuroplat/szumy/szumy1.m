% Skrypt do analizy wynikow z pomiarow szumow. Pomiar za pomoca progqramu w Labview "noise2a", lub "noise2", z katalogu "listopad2003/testy" w "moich dokumentach na Windows. 
%Plik zawiera wyestymowana gestosc widmowa. Fp=20000. N=32768. Normalizacja: FFT, modul kwadrat, razy 4, podzielic przez N, razy fp, razy 1e6. Aby obliczyc sigme, trzeba wrocic do danych przed normalizacja, posumowac CALE widmo (a nie pol), spierwiastkowac, podzielic przez N.

cd /home/pawel/pliki/nauka/neuroplat;
a=importdata('szumy_i_offsety_DAC_ilf20_ihf13_gain24_Vpol050.dat');
%a=importdata('szumy_dummy_DAC_ilf20_ihf13_gain24_Vpol050.dat');
size(a)

fp=20000;
N=32768;
f=[1:16383]/N*fp;

figura1=1;
figura2=2;
figura3=4;

figure(figura1);
for i=1:64
	%figure(figura1);
	%subplot(8,8,i);
	%loglog(f,a(i+1,2:16384)*1e-6);
	%grid on;	
end

figure(figura3);
fontsize=14;

s=a(35,:)*1e-6;
si=sqrt(sum(s)/4*N*fp)/N;
si0=si/900*1e6;
loglog(f,s(1,2:16384)/900/900);
axis([1 10000 1e-16 1e-12])'
h=text(400,1.4e-13,['sigma=' sprintf('%4.2f',si0) 'uV']);
set(h,'FontSize',fontsize);
h=xlabel('frequency[Hz]');
set(h,'FontSize',fontsize);
h=ylabel('power spectrum density [V^{2}/Hz]');
set(h,'FontSize',fontsize);
grid on;

figure(figura2);
plot(a(:,1))
grid on;
