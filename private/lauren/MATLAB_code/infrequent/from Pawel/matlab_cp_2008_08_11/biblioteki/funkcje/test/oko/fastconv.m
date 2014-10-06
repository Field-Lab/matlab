function z=fastconv(a,b,wykladnik);
%Liczenie splotów za pomoc¹ algorytmu FFT. Zwraca wektor o dlugosci
%rownej length(a)-length(b)+1 (jezeli length(a)>length(b)) i s¹
%to dokladnie te probki, dla ktorych nie wystepuja efekty brzegowe
%zwiazane z uzupelnianiem wektorow danych zerami (lub czymkolwiek
%innym). Oznacza to, ¿e na przyklad dla wektorów o identycznych dlugoœciach
%funkcja zwraca tylko jedna probke.
%Maksymalna zajetosc pamieci:
%a) (length(a)+length(b))*2*8 [bytes] lub
%b) max(length(a),lenth(b))*2*8 + 6*(2^wykladnik)*8 [bytes]
%plus ok. 200 bytes.

la=size(a);
if la(1)>1 
   a=a';
end

lb=size(b);
if lb(1)>1 
   b=b';
end

%x - "sygnal"
%h - "filtr" (krotszy niz sygnal)
if (la(2)>lb(2))
   x=a;     
   h=b;
else
   x=b;
   h=a;
end
%whos;
clear a;
clear b;
lx=size(x);
lh=size(h);

z=zeros(1,lx(2)-lh(2)+1);
lz=size(z);

while 2^wykladnik>lz(2)
   wykladnik=wykladnik-1;
end

okno=2^wykladnik;        %ilosc probek dla fft
if okno<lh(2)
   okno=lh(2);
end							 %okno  nie moze byc krotsz niz "flitr"

okno_2=okno-lh(2)+1;      %ilosc prawidlowych probek dla jednego
							    %cyklu obliczen
                         
ilosc=floor(lz(2)/okno_2); %ilosc obiegow petli, z wylaczeniem
								   %ostatniej, obejmujacej probki lokalizoawne
								   %w specjalny sposob
                         
%[fh,r1]=filtr_faza_f(h,okno);
%clear r1;                           
h0=zeros(1,okno);
h0(1,1:lh(2))=h;
fh=fft(h0);
clear h0;
size(fh);

%whos;
for i=1:ilosc
   %i
   poczatek=(i-1)*okno_2;
   e=x(1,(poczatek+1):(poczatek+okno));
   fx=fft(e);
   clear e;
   fx=fx.*fh;
   e0=ifft(fx);
   clear fx;
   e0=real(e0);
   z(1,(poczatek+1):(poczatek+okno_2))=e0(lh(2):okno);
   clear e0;
end
%whos;
%dodatkowy przebieg dla ostatniego okna czasowego
poczatek=lx(2)-okno;
e=x(1,(poczatek+1):(poczatek+okno));
%whos;
fx=fft(e);
%whos;
clear e;
fx=fx.*fh;
%whos;
e0=ifft(fx);
%whos;
clear fx;
e0=real(e0);
z(1,(poczatek+1):(poczatek+okno_2))=e0(lh(2):okno);
%whos;
z=real(z);