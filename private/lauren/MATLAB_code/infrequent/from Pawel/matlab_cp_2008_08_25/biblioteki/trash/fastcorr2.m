function z=fastcorr2(a,b,zasieg,wykladnik);
%z=fastcorr2(a,b,zasieg,wykladnik);
%Liczenie korelacji za pomoca algorytmu FFT. Dane wejsciowe:
%- a,b - wektory danych o IDENTYCZNYCH DLUGOSCIACH;
%- zasieg - tablica dwooch elementow nieujemnych, okreslajaca 
%zadany zasieg liczonych korelacji. Pierwszy element okresla
%opoznienie wektora "b" wzgledem "a", drugi - na odwrot. (Patrz
%dalej;
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


la=size(a);
if la(1)>1 
   a=a';
   la=la';
end

lb=size(b);
if lb(1)>1 
   b=b';
   lb=bl';
end

%x - "sygnal"
%h - "filtr" (krotszy niz sygnal)

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

%la(2)-zasieg(1)-zasieg(2);
ilosc=floor((la(2)-zasieg(1)-zasieg(2))/okno_2);
f=zeros(1,okno);
%h0=zeros(1,okno);

for i=1:ilosc
    x2=zeros(1,okno);
    
    pocz1=(i-1)*okno_2+1;
    koniec1=(i-1)*okno_2+okno;
    pocz2=pocz1+zasieg(1);
    koniec2=koniec1-zasieg(2);
    
    x1=a(1,pocz1:koniec1);
    f1=fft(x1);
    clear x1;
    
    x2(1,1:okno_2)=b(1,koniec2:-1:pocz2);
    f2=fft(x2);
    clear x2;
    
    f1=f1.*f2;
    clear f2;
    e0=ifft(f1);
    %f=f+f1;
    clear f1;
    e0=real(e0);
    %size(e0)
    %plot(e0)
    
    %z0=e0(1,(koniec2-zasieg(1)+1):koniec1);  %bzdura
    %z0=e0(1,(okno-zasieg(1)-zasieg(2)):okno);
    %size(z0)
    z=z+e0(1,(okno-zasieg(1)-zasieg(2)):okno);
end

%x2=zeros(1,okno);

pocz1=koniec1-zasieg(1)-zasieg(2)+1;
koniec1=la(2);
pocz2=pocz1+zasieg(1);
koniec2=koniec1-zasieg(2);

x1=a(1,pocz1:koniec1);
%size(x1);
f1=fft(x1);
clear x1;
    
x2=zeros(1,(koniec1-pocz1+1));
x2(1,1:(koniec2-pocz2+1))=b(1,koniec2:-1:pocz2);
%size(x2);
f2=fft(x2);
clear x2;
    
f1=f1.*f2;
clear f2;
e0=ifft(f1);
clear f1;
e0=real(e0);
se0=size(e0);

z=z+e0(1,(se0(2)-zasieg(1)-zasieg(2)):se0(2));
