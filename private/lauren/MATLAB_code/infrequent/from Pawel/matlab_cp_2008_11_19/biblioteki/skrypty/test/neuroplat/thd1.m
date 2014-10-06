%Pomiar znieksztalcen harmonicznych: sygnal sinusoidalny z generatora, o zmieniajacej sie amplitudie, podawany na wejscie przez dzielnik na plytce (tak jak w "liniowosc1.m"); rejestracja amplitudy sinusa na wej¶ciu oraz znieksztalcen definiowanych jako suma z kwadratem czterech pierwszych harmonicznych (2,3,4,5-krotnosc czest. podtswaowej), na koncu spierwiastkowana. Pomiar dla czestosci sinusa 140Hz i 420Hz. Dla malych amplitud sygnalu widoczny na widmie szum, mimo ze FFT robione z fragmentu danych o dlugosci 5 sekund, zaburza pomiar energii harmonicznych - sta plaska czec nas poczatku. Szum z gruba 1/f, dlatego wieksze wartosci dla malych czest. w przypadku f=140Hz.
%Ustawienia: ilf=20, ihf=13, igain=24, Vpol=0.35

input=[0.1:0.1:3.5]*440e-6;
out_140Hz=[19 39 58 77 96 115 135 154 173 192 211 230 248 267 286 304 323 341 360 378 395 413 430 447 464 482 498 514 529 544 559 573 587 600 613];
thd_140Hz=[0.6 0.6 0.6 0.7 0.8 1.1 1.3 1.7 1.9 2.1 2.5 3.0 3.6 4.2 4.8 5.5 6.1 6.8 7.6 8.5 9.4 10.4 11.4 12.6 13.8 15.3 16.8 18.6 20.6 22.8 25.3 27.9 30.9 34.4 37.9]./out_140Hz*100;
out_420Hz=[22 44 66 88 110 132 154 176 198 219 241 262 284 305 326 347 367 388 408 428 447 467 485 504 522 540 557 573 589 604 618 631 644 656 668];
thd_420Hz=[0.3 0.3 0.3 0.3 0.3 0.4 0.5 0.5 0.7 0.9 1.0 1.3 1.5 1.9 2.3 2.6 3.2 3.9 4.8 5.8 6.9 8.3 10.0 11.9 14.1 16.7 19.5 22.5 26.1 30.0 34.4 39.0 44.2 49.6 55.5]./out_420Hz*100;

figure(1);
clf;
%plot(input,out_140Hz);
h1=gca;
set(h1,'YAxisLocation','right');
hold on;
plot(input*1000,thd_140Hz,'bd-');
xlabel('input P-P [mV]');
ylabel('THD [%]');
h1=gca;
set(h1,'YAxisLocation','left');
title('THD vs input signal level, f=140Hz');
axis([0 1.6 0 10]);
%axis([0 0.0016 0 100]);
grid on;

figure(2);
clf;
%plot(input,out_140Hz);
h1=gca;
set(h1,'YAxisLocation','right');
hold on;
plot(input*1000,thd_420Hz,'bd-');

fontsize=18;

xlabel('input P-P [mV]','FontSize',18);
ylabel('THD [%]','FontSize',18);
%title('THD vs input signal level, f=420Hz');
h1=gca;
set(h1,'YAxisLocation','left');
set(h1,'FontSize',18);
axis([0 1.6 0 10]);
grid on;

dane=[input' thd_420Hz']';
f=fopen('thd.txt','w');
fprintf(f,'%8.6f %8.6f\n',dane);
fclose(f);
q=importdata('thd.txt');
q
size(q)
size(dane)
