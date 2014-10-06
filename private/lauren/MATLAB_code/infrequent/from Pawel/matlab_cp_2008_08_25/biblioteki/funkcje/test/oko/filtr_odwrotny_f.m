function h=filtr_odwrotny_f(filtr,rzad,dynamika);
%Uwaga! Wymagane jest, aby rzad filtru wejsciowego byl nieparzysty.
%dynamika - okresla, ponizej jakiej wartosci wspolczynnika 
%Fouriera (wzgledem wartosci maksymalnej) mozna dany wspolczynnik 
%zamienic na inna wartosc niezerowa (rowna maks/dynamika);
af=size(filtr);

if (af(1)>1)
   filtr=filtr';
end

if ((round(af(2)/2)*2)==af(2))
   error('Rzad filtru wejsciowego musi byc nieparzysty');
end

ar=rzad;

if ((round(rzad/2)*2)==rzad)
   error('Rzad filtru wyjsciowego musi byc nieparzysty');
end

opoznienie_r=floor(af(2)/2);
opoznienie_h=floor(ar/2);

if opoznienie_h<opoznienie_r
   opoznienie_h=opoznienie_r;
   warning('Wydluzona odpowiedz impulsowa filtru wyjsciowego');
end

dlugosc_g=opoznienie_h*2 + 1;
g=zeros(1,dlugosc_g);
g(1,(opoznienie_h+1-opoznienie_r):(opoznienie_h+1+opoznienie_r))=filtr;
ag=size(g);

fg=fft(g);
min_fg=max(max(abs(fg)));
prog_fg=min_fg/dynamika;

for i=1:ag(2)
   if abs(fg(1,i))<prog_fg
      fg(1,i)=prog_fg;
   end
end

ch_cz=fg.^(-1);
%ch_cz=ch_cz.^(-1);
czest=[0:(dlugosc_g-1)]/dlugosc_g;
%ch_cz=ch_cz.*exp(j*2*pi*czest*opoznienie_h);
for i=(opoznienie_h+2):dlugosc_g
   ch_cz(1,i)=conj(ch_cz(1,dlugosc_g+2-i));
end
h=ch_cz;
h=real(ifft(ch_cz));

okno=blackman(rzad)';
h=h.*okno;