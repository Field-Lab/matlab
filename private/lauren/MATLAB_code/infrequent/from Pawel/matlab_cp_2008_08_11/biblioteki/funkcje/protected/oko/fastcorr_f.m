function z=fastcorr_f(read_param,zasieg,wykladnik);
%z=fastcorr_f(read_param,zasieg,wykladnik);
%Estymacja funkcji korelacji skosnej dla dwoch rekordow danych.
%Funkcja samodzielnie odczytuje dane z pliku na podstawie 
%parametrow wejsciowych:
%read_param - struktura okreslajaca dane dla funkcji ->readcnst.
%Pole "channels" zawiera dwa elementy, odpowiadajace dwom kanalom
%z danymi (odpowiada to odpowiednio kanalom "a" oraz "b" dla
%funkcji ->fastcorr);
%- zasieg - tablica dwoch elementow nieujemnych, okreslajaca 
%zadany zasieg liczonych korelacji. Pierwszy element okresla
%opoznienie wektora "b" wzgledem "a", drugi - na odwrot. (Patrz
%dalej);
%- wykladnik - liczba probek na jedno okno czasowe (dla FFT). Op-
%tymalna wartosc tego parametru dla danej dlugosci danych i danego
%zasiegu mozna wyznaczyc za pomoca funkcji... no, na razie jeszcze
%nie ma tej funkcji.
%Wartosci funkcji korelacji dla kazdego wzglednego przesuniecia
%wektorow "a" i "b" sa wyznaczane przy uzyciu tej samej ilosci
%probek. Ilosc ta jest rowna: dlugosc "a" (lub "b") minus zasieg(1)
%minus zasieg(2). W praktyce oznacza to odrzucenie pewnej ilosci 
%probek z wektora "b": poczatkowych "zasieg(1)" probek oraz konco-
%wych "zasieg(2)" probek. 

name=read_param.name;
header=read_param.header;
nrchns=read_param.nrchns;
channels=read_param.channels;
samples=read_param.samples;

la=[1 samples(2)-samples(1)+1];
if la(1)>1 
   a=a';
   la=la';
end

lb=la;
if lb(1)>1 
   b=b';
   lb=bl';
end

if (la(2)~=lb(2))
    error('Dlugosci wektorow nie sa identyczne');
end

while 2^wykladnik>la(2)
   wykladnik=wykladnik-1;
end

okno=2^wykladnik;       %ilosc probek dla fft

if ((zasieg(1)+zasieg(2))>=okno)    
    zasieg=zasieg(1)+zasieg(2)
    okno=okno
    error('Dlugosc okna nie moze byc mniejsza niz calkowity zasieg liczonych korelacji');
end

okno_2=okno-zasieg(1)-zasieg(2);
z=zeros(1,zasieg(1)+zasieg(2)+1);
lz=size(z);

ilosc=floor((la(2)-zasieg(1)-zasieg(2))/okno_2);
f=zeros(1,okno);

for i=1:ilosc
    x2=zeros(1,okno);
    i;
    pocz1=(i-1)*okno_2+samples(1);
    koniec1=pocz1+okno-1;
    pocz2=pocz1+zasieg(1);
    koniec2=koniec1-zasieg(2);
    
    read_param.samples=[pocz1 koniec1];
    read_param.channels=channels(1);
    x1=readcnst(read_param);
    size(x1);
    f1=fft(x1);
    clear x1;
    
    read_param.samples=[pocz2 koniec2];
    read_param.channels=channels(2);
    x2(1,okno_2:-1:1)=readcnst(read_param);
    size(x2);
    f2=fft(x2);
    clear x2;
    
    f1=f1.*f2;
    clear f2;
    f=f+f1;
    clear f1;
end
e0=ifft(f);
clear f;
e0=real(e0);
z=e0(1,(okno-zasieg(1)-zasieg(2)):okno);

pocz1=koniec1-zasieg(1)-zasieg(2)+1;
koniec1=la(2);
pocz2=pocz1+zasieg(1);
koniec2=koniec1-zasieg(2);

read_param.samples=[pocz1 koniec1];
read_param.channels=channels(1);
x1=readcnst(read_param);
size(x1);
f1=fft(x1);
clear x1;
    
x2=zeros(1,(koniec1-pocz1+1));
read_param.samples=[pocz2 koniec2];
read_param.channels=channels(2);
x2(1,(koniec2-pocz2+1):-1:1)=readcnst(read_param);
size(x2);
f2=fft(x2);
clear x2;
    
f1=f1.*f2;
clear f2;
e0=ifft(f1);
clear f1;
e0=real(e0);
se0=size(e0);

z=z+e0(1,(se0(2)-zasieg(1)-zasieg(2)):se0(2));
