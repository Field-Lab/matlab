function [p,Eteor]=maxfit(x,y);
%Funkcja dokonuje fitowania podanego przebiegu (x - argument, y -
%wartosci) funckja liniowa y=ax (nie y=ax+b). Minimalizowany jest blad
%maksymalny, z wykorzystaniem zalozenia, ze przebieg y jest funkcja
%kwadratowa lub ze takie przyblizenie jest akceptowalne (!).
%wyjscia:
%p: tablica postaci [a 0], wprost do wykorzystania z polyval;
%Eteor: wartosc teoretyczna bledu (wzglednego, tzn. podzielonego przez wartosc maksymalna y)
%Oczekuje sie, ze przebieg bledu
%uzyskanego fitu ma identyczne (z dokladnoscia do znaku) wartosci minimalna
%i maksymalna, zblizone do Eteor. Jesli tak nie jest, to przyblizenie
%przebiegu y funkcja kwadratowa jest prawdopodobnie niedokladne.

sy=length(y);

y0=y(sy);
y1=y(round(sy/2));

a2=(2*y0-4*y1)/(x(sy)*x(sy));
a1=(-y0+4*y1)/x(sy);

c=a1+2*a2*x(sy)*(sqrt(2)-1);
p=[c 0];

Eteor=abs((y0-2*y1)*0.35/y0);